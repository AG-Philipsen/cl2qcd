#include "host_kappa.h"

void kappa_karsch (hmc_gaugefield* field, hmc_float & kappa, const hmc_float beta)
{
  //Initialize kappa
  kappa = .0;
  //Compute diagonal spatial components of the energy-momentum-tensor
  hmc_float tdiag_11 [VOL4D];
  hmc_float tdiag_22 [VOL4D];
  hmc_float tdiag_33 [VOL4D];
  //a = 2T_11 - T_22 _ T_33
  hmc_float a [VOL4D];
  //b = 2T_22 - T_11 _ T_33
  hmc_float b [VOL4D];
  //c = 2T_33 - T_22 _ T_11
  hmc_float c [VOL4D];
    
  for (int t=0; t<NTIME; t++){
	for (int n=0; n<VOLSPACE; n++){
	      //Compute required plaquettes
	      hmc_su3matrix temp;
	      
	      local_plaquette(field, & temp, n, t, 1, 0);
	      //faster, to take real first and then trace, but no method
	      hmc_float plaq_10 = trace_su3matrix(&temp).re;
	      local_plaquette(field, & temp, n, t, 2, 0);
	      hmc_float plaq_20 = trace_su3matrix(&temp).re;
	      local_plaquette(field, & temp, n, t, 3, 0);
	      hmc_float plaq_30 = trace_su3matrix(&temp).re;
	      local_plaquette(field, & temp, n, t, 1, 2);
	      hmc_float plaq_12 = trace_su3matrix(&temp).re;
	      local_plaquette(field, & temp, n, t, 1, 3);
	      hmc_float plaq_13 = trace_su3matrix(&temp).re;
	      local_plaquette(field, & temp, n, t, 3, 2);
	      hmc_float plaq_32 = trace_su3matrix(&temp).re;
	      
	      int point = n + VOLSPACE * t;
	      
	      tdiag_11 [point] = plaq_10 + plaq_12 + plaq_13 - plaq_20 - plaq_30 - plaq_32;
	      tdiag_22 [point] = plaq_20 + plaq_12 + plaq_32 - plaq_10 - plaq_30 - plaq_13;
	      tdiag_33 [point] = plaq_30 + plaq_13 + plaq_32 - plaq_10 - plaq_20 - plaq_12;
	
	      a[point] = 2.0*tdiag_11[point] - tdiag_22[point] - tdiag_33[point];
	      b[point] = 2.0*tdiag_22[point] - tdiag_11[point] - tdiag_33[point];
	      c[point] = 2.0*tdiag_33[point] - tdiag_22[point] - tdiag_11[point];
  }}

  hmc_float deltak = 2*PI/NSPACE;
  hmc_float result = 0.0;

  
  for (int x_3=0; x_3<NSPACE; x_3++){
    for (int y_3=0; y_3<x_3; y_3++){
      hmc_float factor = 1.0 - cos(deltak*hmc_float(x_3-y_3));
      for (int x_t=0; x_t<NTIME; x_t++){
	for (int y_t=0; y_t<NTIME; y_t++){
	  for (int x_1=0; x_1<NSPACE; x_1++){
	    for (int y_1=0; y_1<NSPACE; y_1++){
	      for (int x_2=0; x_2<NSPACE; x_2++){
	        for (int y_2=0; y_2<NSPACE; y_2++){
		  
		  int coord_x[NDIM];
		  coord_x[0]= x_t;
		  coord_x[1]= x_1;
		  coord_x[2]= x_2;
		  coord_x[3]= x_3;
		  //new method get_tnspace which gives n_x+VOLSPACE*x_t
		  int n_x = get_nspace (coord_x);
		  int point_x = n_x + VOLSPACE*x_t;
		  int coord_y[NDIM];
		  coord_y[0]= y_t;
		  coord_y[1]= y_1;
		  coord_y[2]= y_2;
		  coord_y[3]= y_3;
		  int n_y = get_nspace (coord_y);
		  int point_y = n_y + VOLSPACE*y_t;

		  result += factor * ( tdiag_11 [point_x] * a[point_y]
		                     + tdiag_22 [point_x] * b[point_y]
		                     + tdiag_33 [point_x] * c[point_y]);
  }}}}}}}}

  //Correlator, 2 by Def, 2 by T_12+T_21  3 by T_12+T_13+T_23 -->/12,
  //Volume for y /VOL4D, /2/pi^2*L_z^2 for derivation, *beta^2 / Nc^2 for T_munu, *2 for y_3<x_3
  
  hmc_float norm = hmc_float (NSPACE* NSPACE) / hmc_float (VOL4D) / hmc_float(NC * NC) /12.0 /PI /PI  * beta * beta;

  kappa = norm * result;
}

void local_Q_plaquette(hmc_3x3matrix * prod, hmc_gaugefield * field, int n, int t, int mu, int nu ){
  hmc_su3matrix tmp;
  int newpos;
  
  //1st plaq
  hmc_su3matrix plaq1;
  //u_mu(x)
  get_su3matrix(&plaq1,field,n,t,mu);
  //u_nu(x+mu)
  if(mu==0) {
    int newt = (t+1)%NTIME;
    get_su3matrix(&tmp,field,n,newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(n,mu),t,nu);
  accumulate_su3matrix_prod(&plaq1,&tmp);
  //adjoint(u_mu(x+nu))
  if(nu==0) {
    int newt = (t+1)%NTIME;
    get_su3matrix(&tmp,field,n,newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(n,nu),t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq1,&tmp);
  //adjoint(u_nu(x))
  get_su3matrix(&tmp,field,n,t,nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq1,&tmp);

  //2nd plaq
  hmc_su3matrix plaq2;
  //U_nu(x)
  get_su3matrix(&plaq2,field,n,t,nu);
  //adj (u_mu(x-mu+nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt =  (t+1)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_neighbor(n,nu),newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(newpos, nu),t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq2,&tmp);
  //adj (u_nu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, nu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq2,&tmp);
  //u_mu(x-mu)
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, mu);
  accumulate_su3matrix_prod(&plaq2,&tmp);

  //3rd plaq
  hmc_su3matrix plaq3;
  //adj (u_mu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, mu);
  adjoin_su3matrix(&tmp);
  copy_su3matrix(&plaq3, &tmp); 
  //adj (u_nu(x-mu-nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, newpos,newt,nu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
  get_su3matrix(&tmp,field,get_lower_neighbor(n,nu),newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(newpos, nu),t,nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq3,&tmp);
  //u_mu(x-mu-nu)
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_lower_neighbor(n,nu),newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(newpos, nu),t,mu);
  accumulate_su3matrix_prod(&plaq3,&tmp);
  //u_nu(x-nu)
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, n,newt, nu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(n, nu),t,nu);
  accumulate_su3matrix_prod(&plaq3,&tmp);

  //4th plaq
  hmc_su3matrix plaq4;
  //adj(u_nu(x-nu))
  adjoin_su3matrix(&tmp);
  copy_su3matrix(&plaq4, &tmp); 
  //u_mu(x-nu)
   if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, nu),t,mu);
  accumulate_su3matrix_prod(&plaq4,&tmp);
  //u_nu(x+mu-nu)
  newpos = get_lower_neighbor(n, nu);
  if (mu==0){
    int newt =  (t+1)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,nu);
  }
  else if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_neighbor(n,mu),newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(newpos, mu),t,nu);
  accumulate_su3matrix_prod(&plaq4,&tmp);
  //adj (u_mu(x))
  get_su3matrix(&tmp,field,n,t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq4,&tmp);

  //Sum up
  su3matrix_to_3x3matrix (prod, &plaq1);
  accumulate_su3matrix_3x3_add(prod, &plaq2);
  accumulate_su3matrix_3x3_add(prod, &plaq3);
  accumulate_su3matrix_3x3_add(prod, &plaq4);
  
}

void testing_Qplak (hmc_complex* plaq, hmc_gaugefield * field)
{
  (*plaq).re=0;
  (*plaq).im=0;
  
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	  for(int nu=0; nu<mu; nu++) {
	    hmc_3x3matrix prod;
	    local_Q_plaquette(&prod, field, n, t, mu, nu );
	    hmc_complex tmpfloat;
	    trace_3x3matrix(&tmpfloat, &prod);
	    (*plaq).re += tmpfloat.re;
	    (*plaq).im += tmpfloat.im;
	    }
	  }
	}
      }
    
  //Divide by 4.0, since Qplaquette consists of 4 plaquettes
  (*plaq).re  *= 2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC) /4.0;
  (*plaq).im  *= 2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC) /4.0;
}


void kappa_clover (hmc_gaugefield* field, hmc_float & kappa, const hmc_float beta){
 
//Initializing result
  kappa = .0;
  
  //Energy-momentum-tensor in clover-discretization
  hmc_float t_12 [VOL4D];
  hmc_float t_13 [VOL4D];
  hmc_float t_23 [VOL4D];
  
  for (int t=0; t<NTIME; t++){
    for (int n=0; n<VOLSPACE; n++){
      //Compute required plaquettes
      hmc_3x3matrix Q_22;
      local_Q_plaquette(&Q_22, field, n, t, 2, 2);
      hmc_3x3matrix Q_10;
      local_Q_plaquette(&Q_10, field, n, t, 1, 0);
      hmc_3x3matrix Q_20;
      local_Q_plaquette(&Q_20, field, n, t, 2, 0);
      hmc_3x3matrix Q_02;
      adjoint_3x3matrix (&Q_02, &Q_20);
// 	      local_Q_plaquette(field, &Q_02, n, t, 0, 2);
      hmc_3x3matrix Q_21;
      local_Q_plaquette(&Q_21,field, n, t, 2, 1);
      hmc_3x3matrix Q_12;
      adjoint_3x3matrix (&Q_12, &Q_21);
// 	      local_Q_plaquette(field, &Q_12, n, t, 1, 2);
      hmc_3x3matrix Q_03;
      local_Q_plaquette(&Q_03, field, n, t, 0, 3);
      hmc_3x3matrix Q_30;
      adjoint_3x3matrix (&Q_30, &Q_03);
      hmc_3x3matrix Q_13;
      local_Q_plaquette( &Q_13, field, n, t, 1, 3);
      hmc_3x3matrix Q_31;
      adjoint_3x3matrix (&Q_31, &Q_13);
      hmc_3x3matrix Q_23;
      local_Q_plaquette( &Q_23, field, n, t, 2, 3);
      hmc_3x3matrix Q_32;
      adjoint_3x3matrix (&Q_32, &Q_23);
// 	      local_Q_plaquette(field, &Q_32, n, t, 3, 2);
      hmc_3x3matrix Q_11;
      local_Q_plaquette( &Q_11, field, n, t, 1, 1);

      int point = n + VOLSPACE * t;
      
      hmc_3x3matrix tmp;
      hmc_complex tmp_cmp;
	      
      //T_12
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_20, &Q_02);
      multiply_3x3matrix (&tmp, &Q_10, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_21, &Q_12);
      multiply_3x3matrix (&tmp, &Q_11, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] += tmp_cmp.re;
      //alpha=2, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_22, &Q_22);
// 	      multiply_3x3matrix (&tmp, &Q_12, &tmp);
// 	      trace_3x3matrix (tmp_cmp, &tmp);
// 	      t_12 [point] += tmp_cmp.re;
      //alpha=3
      subtract_3x3matrix (&tmp, &Q_23, &Q_32);
      multiply_3x3matrix (&tmp, &Q_13, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] += tmp_cmp.re;

      //T_13
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_30, &Q_03);
      multiply_3x3matrix (&tmp, &Q_10, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_31, &Q_13);
      multiply_3x3matrix (&tmp, &Q_11, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] += tmp_cmp.re;
      //alpha=2, vanishes
      subtract_3x3matrix (&tmp, &Q_32, &Q_23);
      multiply_3x3matrix (&tmp, &Q_12, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] += tmp_cmp.re;
      //alpha=3, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_33, &Q_33);
// 	      multiply_3x3matrix (&tmp, &Q_13, &tmp);
// 	      trace_3x3matrix (&tmp_cmp, &tmp);
// 	      t_13 [point] += tmp_cmp.re;

      //T_23
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_30, &Q_03);
      multiply_3x3matrix (&tmp, &Q_20, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_31, &Q_13);
      multiply_3x3matrix (&tmp, &Q_21, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] += tmp_cmp.re;
      //alpha=2, vanishes
      subtract_3x3matrix (&tmp, &Q_32, &Q_23);
      multiply_3x3matrix (&tmp, &Q_22, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] += tmp_cmp.re;
      //alpha=3, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_33, &Q_33);
// 	      multiply_3x3matrix (&tmp, &Q_23, &tmp);
// 	      trace_3x3matrix (&tmp_cmp, &tmp);
// 	      t_23 [point] += tmp_cmp.re;
      }
    }
    
    //Momentum
    const hmc_float deltak = 2.0*PI / NSPACE;
    hmc_float result = 0.0;

  
    for (int x_3=0; x_3<NSPACE; x_3++){
      for (int y_3=0; y_3<x_3; y_3++){
	hmc_float factor = 1.0 - cos(deltak*hmc_float(x_3-y_3));
	for (int x_t=0; x_t<NTIME; x_t++){
	  for (int y_t=0; y_t<NTIME; y_t++){
	    for (int x_1=0; x_1<NSPACE; x_1++){
	      for (int y_1=0; y_1<NSPACE; y_1++){
		for (int x_2=0; x_2<NSPACE; x_2++){
		  for (int y_2=0; y_2<NSPACE; y_2++){
		  
		    int coord_x[NDIM];
		    coord_x[0]= x_t;
		    coord_x[1]= x_1;
		    coord_x[2]= x_2;
		    coord_x[3]= x_3;
		    //new method get_tnspace which gives n_x+VOLSPACE*x_t
		    int n_x = get_nspace (coord_x);
		    int point_x = n_x + VOLSPACE*x_t;
		    int coord_y[NDIM];
		    coord_y[0]= y_t;
		    coord_y[1]= y_1;
		    coord_y[2]= y_2;
		    coord_y[3]= y_3;
		    int n_y = get_nspace (coord_y);
		    int point_y = n_y + VOLSPACE*y_t;

		    //(T_12(x) T_12(y) + T_21(x) T_21(y) + T_13(x) T_13(y)) * factor
		    result += factor* ( t_12[point_x]*t_12[point_y]
				      + t_13[point_x]*t_13[point_y]
				      + t_23[point_x]*t_23[point_y]);
  }}}}}}}}
  //Normalization
  // / (-96)^2=/ 9216, * beta^2, /3 für T_12T_12+T_13T_13+T_23T_23, /VOL4D für Summe y, /2 /pi^2 * L_z^2 für Ableitung * 2 für forschleife
  hmc_float norm = hmc_float (NSPACE* NSPACE) / hmc_float (VOL4D) /PI /PI  * beta * beta /9216. /3.;
    
  kappa = norm * result * beta;

}
