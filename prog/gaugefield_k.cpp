#include "gaugefield_k.h"

hmc_error Gaugefield_k::init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer){
	
	hmc_error err;
  
	hmc_gaugefield * gf_tmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	err = set_gf (gf_tmp);
	
	err = set_parameters (input_parameters);

	init_gaugefield(timer);

	err = set_num_ocl_devices (numdevs);

	if(numdevs != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0){
		Opencl_k * devices_tmp = new Opencl_k[get_num_ocl_devices()];
		err = set_devices (devices_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		cout << "init device #" << n << endl;
		(get_devices())[n].init(devicetypes[n], local_work_size, global_work_size, timer, get_parameters ());
		
	}

	return HMC_SUCCESS;
}

hmc_float Gaugefield_k::get_kappa_karsch (){
	return kappa_karsch_val;
}
	
hmc_float Gaugefield_k::get_kappa_clover (){
	return kappa_clover_val;
}

hmc_error Gaugefield_k::set_kappa_karsch (hmc_float in){
	kappa_karsch_val = in;
	return HMC_SUCCESS;
}

hmc_error Gaugefield_k::set_kappa_clover (hmc_float in){
	kappa_clover_val = in;
	return HMC_SUCCESS;
  
}

hmc_error Gaugefield_k::kappa_karsch_gpu (const size_t local_work_size, const size_t global_work_size, usetimer* timer_karsch){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = get_devices_k()[0].run_kappa_karsch_gpu(local_work_size, global_work_size, timer_karsch);
	return err;
}

hmc_error Gaugefield_k::kappa_clover_gpu (const size_t local_work_size, const size_t global_work_size, usetimer* timer_clover){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = get_devices_k()[0].run_kappa_clover_gpu(local_work_size, global_work_size, timer_clover);
	return err;
}

Opencl * Gaugefield_k::get_devices_k (){
  return  devices_k;
}

hmc_error Gaugefield_k::kappa_karsch ()
{
  //Initialize kappa
  kappa_karsch_val = 0.0;
  //Initializing beta
  const hmc_float beta = get_parameters()->get_beta();
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
	      
	      local_plaquette(get_gf(), & temp, n, t, 1, 0);
	      //faster, to take real first and then trace, but no method
	      hmc_float plaq_10 = trace_su3matrix(&temp).re;
	      local_plaquette(get_gf(), & temp, n, t, 2, 0);
	      hmc_float plaq_20 = trace_su3matrix(&temp).re;
	      local_plaquette(get_gf(), & temp, n, t, 3, 0);
	      hmc_float plaq_30 = trace_su3matrix(&temp).re;
	      local_plaquette(get_gf(), & temp, n, t, 1, 2);
	      hmc_float plaq_12 = trace_su3matrix(&temp).re;
	      local_plaquette(get_gf(), & temp, n, t, 1, 3);
	      hmc_float plaq_13 = trace_su3matrix(&temp).re;
	      local_plaquette(get_gf(), & temp, n, t, 3, 2);
	      hmc_float plaq_32 = trace_su3matrix(&temp).re;

	      int point = n + VOLSPACE * t;
	      
	      tdiag_11 [point] = plaq_10 + plaq_12 + plaq_13 - plaq_20 - plaq_30 - plaq_32;
	      tdiag_22 [point] = plaq_20 + plaq_12 + plaq_32 - plaq_10 - plaq_30 - plaq_13;
	      tdiag_33 [point] = plaq_30 + plaq_13 + plaq_32 - plaq_10 - plaq_20 - plaq_12;

	      a[point] = 2.0*tdiag_11[point] - tdiag_22[point] - tdiag_33[point];
	      b[point] = 2.0*tdiag_22[point] - tdiag_11[point] - tdiag_33[point];
	      c[point] = 2.0*tdiag_33[point] - tdiag_22[point] - tdiag_11[point];
  }}

  hmc_float deltak = 2*PI/ hmc_float(NSPACE);
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

  kappa_karsch_val = norm * result;
  
  return HMC_SUCCESS;
}

hmc_error Gaugefield_k::kappa_clover (){
  //Initializing result
  kappa_clover_val = .0;
  //Initializing beta
  const hmc_float beta = get_parameters()->get_beta();
  //Energy-momentum-tensor in clover-discretization
  hmc_float t_12 [VOL4D];
  hmc_float t_13 [VOL4D];
  hmc_float t_23 [VOL4D];
  
  for (int t=0; t<NTIME; t++){
    for (int n=0; n<VOLSPACE; n++){
      //Compute required plaquettes
      hmc_3x3matrix Q_22;
      local_Q_plaquette(&Q_22, get_gf(), n, t, 2, 2);
      hmc_3x3matrix Q_10;
      local_Q_plaquette(&Q_10, get_gf(), n, t, 1, 0);
      hmc_3x3matrix Q_20;
      local_Q_plaquette(&Q_20, get_gf(), n, t, 2, 0);
      hmc_3x3matrix Q_02;
      adjoint_3x3matrix (&Q_02, &Q_20);
      hmc_3x3matrix Q_21;
      local_Q_plaquette(&Q_21,get_gf(), n, t, 2, 1);
      hmc_3x3matrix Q_12;
      adjoint_3x3matrix (&Q_12, &Q_21);
      hmc_3x3matrix Q_03;
      local_Q_plaquette(&Q_03, get_gf(), n, t, 0, 3);
      hmc_3x3matrix Q_30;
      adjoint_3x3matrix (&Q_30, &Q_03);
      hmc_3x3matrix Q_13;
      local_Q_plaquette( &Q_13, get_gf(), n, t, 1, 3);
      hmc_3x3matrix Q_31;
      adjoint_3x3matrix (&Q_31, &Q_13);
      hmc_3x3matrix Q_23;
      local_Q_plaquette( &Q_23, get_gf(), n, t, 2, 3);
      hmc_3x3matrix Q_32;
      adjoint_3x3matrix (&Q_32, &Q_23);
      hmc_3x3matrix Q_11;
      local_Q_plaquette( &Q_11, get_gf(), n, t, 1, 1);

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
      //alpha=2
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
      //alpha=2
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
    const hmc_float deltak = 2.0*PI / hmc_float (NSPACE);
    hmc_float result = 0.0;

  
    for (int x_3=0; x_3<NSPACE; x_3++){
      for (int y_3=0; y_3<x_3; y_3++){
	hmc_float factor = 1.0 - cos(deltak * hmc_float(x_3-y_3));
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
  // 1/3 for averaging T_ij, 1/V/Nt for averaging y, L^2/2/pi^2 for derivation, (-1/64)^2 for Clover and T_munu^2, beta^2/Nc^2 for T_munu^2
  // *2 for temp + conj (temp) *2 for for-loop
  // = beta^2 * L^2/ (55296 * V * Nt * pi^2)
  hmc_float norm = hmc_float (NSPACE* NSPACE) / hmc_float (VOL4D) /PI /PI * beta * beta / 55296. ;
    
  kappa_clover_val = norm * result;
  
  return HMC_SUCCESS;
}