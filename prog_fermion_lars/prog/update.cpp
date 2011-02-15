#include "update.h"
#include "testing.h"
// Construct new SU2 matrix using improved alg by Kennedy Pendleton
void SU2Update(hmc_float dst [su2_entries], const hmc_float alpha)
{
  hmc_float delta;
  hmc_float a0 ;
  hmc_float eta ; int cter = 0;
  do
  {
    delta = -log(myran.doub())/alpha*pow(cos(2. * PI * myran.doub()), 2.) -log(myran.doub())/alpha;
    a0 = 1.-delta;
    eta = myran.doub();
    cter ++;
  }while ( (1.-0.5*delta) < eta*eta); 
  hmc_float phi = 2.*PI*myran.doub();
  hmc_float theta = asin(2.*myran.doub() - 1.);
  dst[0] = a0;
  dst[1] = sqrt(1.-a0 * a0)*cos(theta) * cos(phi);
  dst[2] = sqrt(1.-a0 * a0)*cos(theta) * sin(phi);
  dst[3] = sqrt(1.-a0 * a0)*sin(theta);
}

// void calc_staple(hmc_gaugefield * field, hmc_su3matrix * dest, const int pos, const int t, const int mu_in){
void calc_staple(hmc_gaugefield * field, hmc_staplematrix * dest, const int pos, const int t, const int mu_in){
  hmc_su3matrix prod, prod2, tmp;
  hmc_staplematrix dummy;
  int nu, newpos, newt;
  
  zero_staplematrix(&dummy);
  //iterate through the three directions other than mu
  for(int i = 1; i<NDIM; i++){
    zero_su3matrix(&prod);
    nu = (mu_in + i)%NDIM;
    //first staple
    //u_nu(x+mu)
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(&tmp,field,pos,newt,nu);
    } else {
      get_su3matrix(&tmp,field,get_neighbor(pos,mu_in),t,nu);
    }
    copy_su3matrix(&prod, &tmp);
    //adjoint(u_mu(x+nu))
    if(nu==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(&tmp,field,pos,newt,mu_in);
    } else {
      get_su3matrix(&tmp,field,get_neighbor(pos,nu),t,mu_in);
    }
    adjoin_su3matrix(&tmp);
    accumulate_su3matrix_prod(&prod,&tmp);
    //adjoint(u_nu(x))
    get_su3matrix(&tmp,field,pos,t,nu);
    adjoin_su3matrix(&tmp);
    accumulate_su3matrix_prod(&prod,&tmp);  
   
    //second staple
    //adjoint (u_nu(x+mu-nu))
    //newpos is "pos-nu" (spatial)
    newpos = get_lower_neighbor(pos, nu);
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(&tmp,field,newpos,newt,nu);
    } else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(&tmp,field,get_neighbor(pos,mu_in),newt,nu);
    }
    else{
	get_su3matrix(&tmp,field,get_neighbor(newpos,mu_in),t,nu);
    }
    adjoin_su3matrix(&tmp);
    copy_su3matrix(&prod2, &tmp);
    //adjoint(u_mu(x-nu))
    if(mu_in==0) {
      get_su3matrix(&tmp,field,newpos,t,mu_in);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(&tmp,field,pos,newt,mu_in);
    }
    else{
	get_su3matrix(&tmp,field,newpos,t,mu_in);
    }
    adjoin_su3matrix(&tmp);
    accumulate_su3matrix_prod(&prod2,&tmp);
    //adjoint(u_nu(x-nu))
    if(mu_in==0) {
      get_su3matrix(&tmp,field,newpos,t,nu);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(&tmp,field,pos,newt,nu);
    }
    else{
	get_su3matrix(&tmp,field,newpos,t,nu);
    }
    accumulate_su3matrix_prod(&prod2,&tmp); 

    accumulate_su3matrices_add(&dummy, &prod);
    accumulate_su3matrices_add(&dummy, &prod2);
  }
  copy_staplematrix(dest, &dummy);
}

void heatbath_update (hmc_gaugefield * gaugefield, const hmc_float beta){
  int pos, mu, t; 
  //iterate through the sites
  for(t=0; t<NTIME; t++){
    for(pos=0;pos<VOLSPACE;pos++){
      for(mu = 0; mu<NDIM; mu++){
      hmc_su3matrix U;
      hmc_staplematrix W,  staplematrix;
      
      hmc_complex w [su2_entries];
      hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
      int order[3]; 
      random_1_2_3(order);
      //old U
      get_su3matrix(&U, gaugefield, pos, t, mu);
      
      calc_staple(gaugefield, &staplematrix, pos, t, mu);
      for(int i=0; i<3; i++)
      {
        //Produkt aus U und staple, saved in W
        multiply_staplematrix(&W, &U, &staplematrix);
        //Reduktion auf SU(2) submatrix
        reduction (w, W, order[i]);
        //Darstellung von W in Paulibasis
        w_pauli[0] = 0.5*(w[0].re + w[3].re);
        w_pauli[1] = 0.5*(w[1].im + w[2].im);
        w_pauli[2] = 0.5*(w[1].re - w[2].re);
        w_pauli[3] = 0.5*(w[0].im - w[3].im);

        //Berechnung der Normierung k, entspricht Determinante in normaler Basis
        k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
        w_pauli[0] = hmc_float(1)/k * w_pauli[0];
        w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
        w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
        w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
        //beta' = 2Beta/Nc*k
        hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
      
        //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
        SU2Update(r_pauli, beta_neu);
        //Multipliziere neuen Link mit w, alles in Paulibasis
	hmc_float su2_tmp[su2_entries];
        su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
        su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
        su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
        su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
        r_pauli[0] = su2_tmp[0];      
        r_pauli[1] = su2_tmp[1];
        r_pauli[2] = su2_tmp[2];
        r_pauli[3] = su2_tmp[3];

      
        //go back to a su2 matrix in standard basis
        w[0].re = r_pauli[0];
        w[0].im = r_pauli[3];
        w[1].re = r_pauli[2];
        w[1].im = r_pauli[1];
        w[2].re = -r_pauli[2];
        w[2].im = r_pauli[1];
        w[3].re = r_pauli[0];
        w[3].im = -r_pauli[3];
  
        //extend to SU3
        hmc_su3matrix extW, tmp;
        extend (&extW, order[i], w);
	
	multiply_su3matrices(&tmp, &extW, &U);

	copy_su3matrix(&U, &tmp);
      }
      put_su3matrix(gaugefield, &U, pos, t, mu);
  }}}
}

void heatbath_overrelax (hmc_gaugefield * gaugefield, const hmc_float beta){
  int pos, mu, t;
  //iterate through the sites
  for(t=0; t<NTIME; t++){
    for(pos=0;pos<VOLSPACE;pos++){
      for(mu = 0; mu<NDIM; mu++){
      hmc_su3matrix U;
      hmc_staplematrix W,  staplematrix;
      hmc_complex w [su2_entries];
      hmc_float w_pauli[su2_entries], k;
      int order[3]; 
      random_1_2_3(order);
      
      //old U
      get_su3matrix(&U, gaugefield, pos, t, mu);
  
      calc_staple(gaugefield, &staplematrix, pos, t, mu);
    
      for(int i=0; i<3; i++)
      {
        //Produkt aus U und staple, saved in W
        multiply_staplematrix(&W, &U, &staplematrix);
	
	//Reduktion auf SU(2) submatrix
        reduction (w, W, order[i]);
	
        //Darstellung von W in Paulibasis
        w_pauli[0] = 0.5*(w[0].re + w[3].re);
        w_pauli[1] = 0.5*(w[1].im + w[2].im);
        w_pauli[2] = 0.5*(w[1].re - w[2].re);
        w_pauli[3] = 0.5*(w[0].im - w[3].im);
     
        //Berechnung der Normierung k, entspricht Determinante in normaler Basis
        k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
        hmc_float a[su2_entries];
        a[0] = hmc_float(1)/k * w_pauli[0];
        a[1] = hmc_float(-1)/k * w_pauli[1];
        a[2] = hmc_float(-1)/k * w_pauli[2];
        a[3] = hmc_float(-1)/k * w_pauli[3];
        //Square a and save in w
        w_pauli[0] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
        w_pauli[1] = 2.*a[0]*a[1];
        w_pauli[2] = 2.*a[0]*a[2];
        w_pauli[3] = 2.*a[0]*a[3];
      
        //go back to a su2 matrix in standard basis
        w[0].re = w_pauli[0];
        w[0].im = w_pauli[3];
        w[1].re = w_pauli[2];
        w[1].im = w_pauli[1];
        w[2].re = -w_pauli[2];
        w[2].im = w_pauli[1];
        w[3].re = w_pauli[0];
        w[3].im = -w_pauli[3];
	
        //extend to SU3
        hmc_su3matrix extW, tmp;
        extend (&extW, order[i], w);
	
	multiply_su3matrices(&tmp, &extW, &U);
	copy_su3matrix(&U, &tmp);
      }
      put_su3matrix(gaugefield, &U, pos, t, mu);
  }}}
}

void heatbath_update_checkerboard (hmc_gaugefield * gaugefield, const hmc_float beta){
  int pos, mu, t, i, x, y, z;
  //iterate through the sites
  for (mu = 0; mu < NDIM; mu ++){
    //update even sites
    #pragma omp for private (pos, i, x, y, z, t)
    for (t = 0; t<NTIME; t++){
      for (x = 0; x<NSPACE; x++){
        for (y = 0; y<NSPACE/2; y++){
          for (z = 0; z< NSPACE; z++){ 
	
	  //!! the formula should be:
	  //!!  mu + NDIM*( int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x ) + NDIM*VOLSPACE*h;
	  //!!from which one gets:
	  pos = int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x;
	
	  hmc_su3matrix U;
          hmc_staplematrix W,  staplematrix;
      
          hmc_complex w [su2_entries];
          hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
          int order[3]; 
          random_1_2_3(order);
          //old U
          get_su3matrix(&U, gaugefield, pos, t, mu);
          
          calc_staple(gaugefield, &staplematrix, pos, t, mu);
          for(i=0; i<3; i++)
          {
            //Produkt aus U und staple, saved in W
            multiply_staplematrix(&W, &U, &staplematrix);
            //Reduktion auf SU(2) submatrix
            reduction (w, W, order[i]);
            //Darstellung von W in Paulibasis
            w_pauli[0] = 0.5*(w[0].re + w[3].re);
            w_pauli[1] = 0.5*(w[1].im + w[2].im);
            w_pauli[2] = 0.5*(w[1].re - w[2].re);
            w_pauli[3] = 0.5*(w[0].im - w[3].im);

            //Berechnung der Normierung k, entspricht Determinante in normaler Basis
            k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
            w_pauli[0] = hmc_float(1)/k * w_pauli[0];
            w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
            w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
            w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
            //beta' = 2Beta/Nc*k
            hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
      
            //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
            SU2Update(r_pauli, beta_neu);
            //Multipliziere neuen Link mit w, alles in Paulibasis
            hmc_float su2_tmp[su2_entries];
            su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
            su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
            su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
            su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
            r_pauli[0] = su2_tmp[0];      
            r_pauli[1] = su2_tmp[1];
            r_pauli[2] = su2_tmp[2];
            r_pauli[3] = su2_tmp[3];

            //go back to a su2 matrix in standard basis
            w[0].re = r_pauli[0];
            w[0].im = r_pauli[3];
            w[1].re = r_pauli[2];
            w[1].im = r_pauli[1];
            w[2].re = -r_pauli[2];
            w[2].im = r_pauli[1];
            w[3].re = r_pauli[0];
            w[3].im = -r_pauli[3];
  
            //extend to SU3
            hmc_su3matrix extW, tmp;
            extend (&extW, order[i], w);
	
	    multiply_su3matrices(&tmp, &extW, &U);

	    copy_su3matrix(&U, &tmp);
          }
          put_su3matrix(gaugefield, &U, pos, t, mu); 
    }}}}
    //make sure all even sites are done
    #pragma omp barrier

    //update odd sites
    #pragma omp for private (pos, i, x, y, z, t)
    for (t = 0; t<NTIME; t++){
      for (x = 0; x<NSPACE; x++){
        for (y = 0; y<NSPACE/2; y++){
          for (z = 0; z< NSPACE; z++){ 
	
	  //!! the formula should be:
	  //!!  mu + NDIM*( int((x+t+1)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x ) + NDIM*VOLSPACE*h;
	  //!!from which one gets:
	  pos = int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x;

          hmc_su3matrix U;
          hmc_staplematrix W,  staplematrix;
      
          hmc_complex w [su2_entries];
          hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
          int order[3]; 
          random_1_2_3(order);
          //old U
          get_su3matrix(&U, gaugefield, pos, t, mu);
          
          calc_staple(gaugefield, &staplematrix, pos, t, mu);
          for(i=0; i<3; i++)
          {
            //Produkt aus U und staple, saved in W
            multiply_staplematrix(&W, &U, &staplematrix);
            //Reduktion auf SU(2) submatrix
            reduction (w, W, order[i]);
            //Darstellung von W in Paulibasis
            w_pauli[0] = 0.5*(w[0].re + w[3].re);
            w_pauli[1] = 0.5*(w[1].im + w[2].im);
            w_pauli[2] = 0.5*(w[1].re - w[2].re);
            w_pauli[3] = 0.5*(w[0].im - w[3].im);

            //Berechnung der Normierung k, entspricht Determinante in normaler Basis
            k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
            w_pauli[0] = hmc_float(1)/k * w_pauli[0];
            w_pauli[1] = hmc_float(-1)/k * w_pauli[1];
            w_pauli[2] = hmc_float(-1)/k * w_pauli[2];
            w_pauli[3] = hmc_float(-1)/k * w_pauli[3];
	
            //beta' = 2Beta/Nc*k
            hmc_float beta_neu =  2.*beta / hmc_float(NC)*k;
      
            //Neuer Link in Paulibasis mittels Kennedy-Pendleton-Algorithmus
            SU2Update(r_pauli, beta_neu);
            //Multipliziere neuen Link mit w, alles in Paulibasis
            hmc_float su2_tmp[su2_entries];
            su2_tmp[0] = r_pauli[0]*w_pauli[0] - r_pauli[1]*w_pauli[1] - r_pauli[2]*w_pauli[2] - r_pauli[3]*w_pauli[3] ;
            su2_tmp[1] = w_pauli[0]*r_pauli[1] + w_pauli[1]*r_pauli[0] - r_pauli[2]*w_pauli[3] + r_pauli[3]*w_pauli[2] ;
            su2_tmp[2] = w_pauli[0]*r_pauli[2] + w_pauli[2]*r_pauli[0] - r_pauli[3]*w_pauli[1] + r_pauli[1]*w_pauli[3] ;
            su2_tmp[3] = w_pauli[0]*r_pauli[3] + w_pauli[3]*r_pauli[0] - r_pauli[1]*w_pauli[2] + r_pauli[2]*w_pauli[1] ;
            r_pauli[0] = su2_tmp[0];      
            r_pauli[1] = su2_tmp[1];
            r_pauli[2] = su2_tmp[2];
            r_pauli[3] = su2_tmp[3];

            //go back to a su2 matrix in standard basis
            w[0].re = r_pauli[0];
            w[0].im = r_pauli[3];
            w[1].re = r_pauli[2];
            w[1].im = r_pauli[1];
            w[2].re = -r_pauli[2];
            w[2].im = r_pauli[1];
            w[3].re = r_pauli[0];
            w[3].im = -r_pauli[3];
  
            //extend to SU3
            hmc_su3matrix extW, tmp;
            extend (&extW, order[i], w);
	
	    multiply_su3matrices(&tmp, &extW, &U);

	    copy_su3matrix(&U, &tmp);
          }
          put_su3matrix(gaugefield, &U, pos, t, mu); 
    }}}}
  }
}

void heatbath_overrelax_checkerboard (hmc_gaugefield * gaugefield, const hmc_float beta){
  int pos, i, mu, t, x, y, z;
  //iterate through the sites
  for (mu = 0; mu < NDIM; mu ++){
    //update even sites
    #pragma omp for private (pos, i, x, y, z, t)
    for (t = 0; t<NTIME; t++){
      for (x = 0; x<NSPACE; x++){
        for (y = 0; y<NSPACE/2; y++){
          for (z = 0; z< NSPACE; z++){ 
	
	  //!! the formula should be:
	  //!!  mu + NDIM*( int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x ) + NDIM*VOLSPACE*h;
	  //!!from which one gets:
	  pos = int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x;

          hmc_su3matrix U;
          hmc_staplematrix W,  staplematrix;
          hmc_complex w [su2_entries];
          hmc_float w_pauli[su2_entries], k;
          int order[3]; 
          random_1_2_3(order);
      
          //old U
          get_su3matrix(&U, gaugefield, pos, t, mu);
  
          calc_staple(gaugefield, &staplematrix, pos, t, mu);
    
          for(i=0; i<3; i++)
          {
            //Produkt aus U und staple, saved in W
            multiply_staplematrix(&W, &U, &staplematrix);
	
	    //Reduktion auf SU(2) submatrix
            reduction (w, W, order[i]);
	
            //Darstellung von W in Paulibasis
            w_pauli[0] = 0.5*(w[0].re + w[3].re);
            w_pauli[1] = 0.5*(w[1].im + w[2].im);
            w_pauli[2] = 0.5*(w[1].re - w[2].re);
            w_pauli[3] = 0.5*(w[0].im - w[3].im);
     
            //Berechnung der Normierung k, entspricht Determinante in normaler Basis
            k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
            hmc_float a[su2_entries];
            a[0] = hmc_float(1)/k * w_pauli[0];
            a[1] = hmc_float(-1)/k * w_pauli[1];
            a[2] = hmc_float(-1)/k * w_pauli[2];
            a[3] = hmc_float(-1)/k * w_pauli[3];
            //Square a and save in w
            w_pauli[0] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
            w_pauli[1] = 2.*a[0]*a[1];
            w_pauli[2] = 2.*a[0]*a[2];
            w_pauli[3] = 2.*a[0]*a[3];
      
            //go back to a su2 matrix in standard basis
            w[0].re = w_pauli[0];
            w[0].im = w_pauli[3];
            w[1].re = w_pauli[2];
            w[1].im = w_pauli[1];
            w[2].re = -w_pauli[2];
            w[2].im = w_pauli[1];
            w[3].re = w_pauli[0];
            w[3].im = -w_pauli[3];
	
            //extend to SU3
            hmc_su3matrix extW, tmp;
            extend (&extW, order[i], w);
	
	    multiply_su3matrices(&tmp, &extW, &U);
	    copy_su3matrix(&U, &tmp);
          }
          put_su3matrix(gaugefield, &U, pos, t, mu);
    }}}}
    //make sure all even sites are done
    #pragma omp barrier

    //update odd sites
    #pragma omp for private (pos, i, x, y, z, t)
    for (t = 0; t<NTIME; t++){
      for (x = 0; x<NSPACE; x++){
        for (y = 0; y<NSPACE/2; y++){
          for (z = 0; z< NSPACE; z++){ 
	
	  //!! the formula should be:
	  //!!  mu + NDIM*( int((x+t+1)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x ) + NDIM*VOLSPACE*h;
	  //!!from which one gets:
	  pos = int((x+t)%2)*(1 + 2*z - int (2*z/NSPACE)) + int((x+1+t)%2)*(2*z + int (2*z/NSPACE)) + 2*NSPACE*y + NSPACE*NSPACE*x;

           hmc_su3matrix U;
           hmc_staplematrix W,  staplematrix;
           hmc_complex w [su2_entries];
           hmc_float w_pauli[su2_entries], k;
           int order[3]; 
           random_1_2_3(order);
      
           //old U
           get_su3matrix(&U, gaugefield, pos, t, mu);
   
           calc_staple(gaugefield, &staplematrix, pos, t, mu);
    
           for(i=0; i<3; i++)
           {
             //Produkt aus U und staple, saved in W
             multiply_staplematrix(&W, &U, &staplematrix);
	 
  	     //Reduktion auf SU(2) submatrix
             reduction (w, W, order[i]);
	
             //Darstellung von W in Paulibasis
             w_pauli[0] = 0.5*(w[0].re + w[3].re);
             w_pauli[1] = 0.5*(w[1].im + w[2].im);
             w_pauli[2] = 0.5*(w[1].re - w[2].re);
             w_pauli[3] = 0.5*(w[0].im - w[3].im);
     
             //Berechnung der Normierung k, entspricht Determinante in normaler Basis
             k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );
             hmc_float a[su2_entries];
             a[0] = hmc_float(1)/k * w_pauli[0];
             a[1] = hmc_float(-1)/k * w_pauli[1];
             a[2] = hmc_float(-1)/k * w_pauli[2];
             a[3] = hmc_float(-1)/k * w_pauli[3];
             //Square a and save in w
             w_pauli[0] = a[0]*a[0] - a[1]*a[1] - a[2]*a[2] - a[3]*a[3];
             w_pauli[1] = 2.*a[0]*a[1];
             w_pauli[2] = 2.*a[0]*a[2];
             w_pauli[3] = 2.*a[0]*a[3];
       
             //go back to a su2 matrix in standard basis
             w[0].re = w_pauli[0];
             w[0].im = w_pauli[3];
             w[1].re = w_pauli[2];
             w[1].im = w_pauli[1];
             w[2].re = -w_pauli[2];
             w[2].im = w_pauli[1];
             w[3].re = w_pauli[0];
             w[3].im = -w_pauli[3];
	
             //extend to SU3
             hmc_su3matrix extW, tmp;
             extend (&extW, order[i], w);
	
 	     multiply_su3matrices(&tmp, &extW, &U);
 	     copy_su3matrix(&U, &tmp);
          }
          put_su3matrix(gaugefield, &U, pos, t, mu);
    }}}}
  }
}