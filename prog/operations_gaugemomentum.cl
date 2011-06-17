/**
 * @file operations on algebraelements
 */

ae set_zero_ae(){
	ae tmp;
	tmp.e0 = 0.;
	tmp.e1 = 0.;
	tmp.e2 = 0.;
	tmp.e3 = 0.;
	tmp.e4 = 0.;
	tmp.e5 = 0.;
	tmp.e6 = 0.;
	tmp.e7 = 0.;
	return tmp;
}

hmc_float ae_squarenorm(ae in){
	hmc_float result = 
		(in).e0*(in).e0 +
		(in).e1*(in).e1 +
		(in).e2*(in).e2 +
		(in).e3*(in).e3 +
		(in).e4*(in).e4 +
		(in).e5*(in).e5 +
		(in).e6*(in).e6 +
		(in).e7*(in).e7;
	return result;
}

ae acc_factor_times_algebraelement(ae in, hmc_float factor, ae force_in){
	ae tmp;
	tmp.e0 = (in).e0+factor*(force_in).e0; 
	tmp.e1 = (in).e1+factor*(force_in).e1;
	tmp.e2 = (in).e2+factor*(force_in).e2;
	tmp.e3 = (in).e3+factor*(force_in).e3;
	tmp.e4 = (in).e4+factor*(force_in).e4;
	tmp.e5 = (in).e5+factor*(force_in).e5;
	tmp.e6 = (in).e6+factor*(force_in).e6;
	tmp.e7 = (in).e7+factor*(force_in).e7;
	return tmp;
}

//calculates the trace of i times generator times 3x3-matrix and stores this in a su3-algebraelement
//now using structs
ae tr_lambda_u(Matrix3x3 in){
	ae tmp;
	tmp.e0 = ( -in.e10.im - in.e01.im);
	tmp.e1 = (+in.e10.re-in.e01.re);
	tmp.e2 = (-in.e00.im+in.e11.im);
	tmp.e3 = (-in.e20.im-in.e02.im);
	tmp.e4 = (+in.e20.re-in.e02.re);
	tmp.e5 = (-in.e21.im-in.e12.im);
	tmp.e6 = (+in.e21.re-in.e12.re);
	tmp.e7 = (-in.e00.im-in.e11.im + 2.0*in.e22.im)*0.577350269189625;
	return tmp;
}

/**
 * @file operations used by gaugemomentum
 */


/** @todo memcpy ... */
/** @todo get rid of the workaround.. */
// hmc_error copy_gaugemomenta(hmc_gauge_momentum * source, hmc_gauge_momentum * dest){
// 	// copies source to destination within cpu memory, layer for momentum array
// 	return hmc_floatcopy((hmc_float *)source, (hmc_float *)dest, GAUGEMOMENTASIZE); // SL: not tested
// }

//deprecated, this can be done by in = out in the code...
/*
hmc_algebraelement copy_gaugemomenta(hmc_algebraelement2 * source, hmc_algebraelement2 * dest){
	hmc_algebraelement2 tmp;
	for(int i = 0; i<GAUGEMOMENTASIZE; i++){	
		(dest[i]).e0 = (source[i]).e0;
		(dest[i]).e1 = (source[i]).e1;
		(dest[i]).e2 = (source[i]).e2;
		(dest[i]).e3 = (source[i]).e3;
		(dest[i]).e4 = (source[i]).e4;
		(dest[i]).e5 = (source[i]).e5;
		(dest[i]).e6 = (source[i]).e6;
		(dest[i]).e7 = (source[i]).e7;
	}
	return HMC_SUCCESS;
}
*/
	
/** @todo add args for reduction... */
__kernel void gaugemomenta_squarenorm(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		hmc_float result = 0.;
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			result += ae_squarenorm(in[i]);
		}
		
		/** @todo add reduction.. */
		
	}
}

/** @todo memset... */
__kernel void set_zero_gaugemomentum(__global ae * in){
	int id = get_global_id(0);
	if(id == 0){
		for(int i = 0; i<GAUGEMOMENTASIZE; i++){
			in[i] = set_zero_ae();
		}
	}
}

/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global ae * out, __global hmc_ocl_ran * rnd){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_complex tmp;
	
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {	
		//CP: there are 8 elements in ae
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e0 = tmp.re;
		out[id_tmp].e1 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e2 = tmp.re;
		out[id_tmp].e3 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e4 = tmp.re;
		out[id_tmp].e5 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e6 = tmp.re;
		out[id_tmp].e7 = tmp.im;
	}
}



void update_gaugemomentum(ae in, hmc_float factor, int global_link_pos, __global ae * out){
	ae tmp = out[global_link_pos];
	acc_factor_times_algebraelement(tmp, factor, in);
	out[global_link_pos] = tmp;
}


