// Description of variables of this kernel:
//  - out: The output staggered field: gaussian complex su3vec (site by site)
//  - rngStates: The state of the generator (rngStateStorageType is useful to produce
//               random number --> see random.cl)

__kernel void set_gaussian_spinorfield_stagg_eoprec(__global staggeredStorageType * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	
	// complex to store the gaussian numbers drawn
	hmc_complex tmp;
	// sigma has to be 0.5 here (for the explanation of why, see the documentation
	// of the function set_gaussian_spinorfield_device in spinors_staggered.hpp)
	hmc_float sigma = 0.5;
	
	su3vec out_tmp;
	
	//if one wants to compare rnd numbers as from a single threaded program
#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	for(int id_local = id; id_local < EOPREC_SPINORFIELDSIZE_LOCAL; id_local += global_size) {
		site_idx id_mem = get_eo_site_idx_from_st_idx(get_even_st_idx_local(id_local));
		
		//There are NC=3 complex elements in the su3vec
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.re = tmp.re;
		out_tmp.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.re = tmp.re;
		out_tmp.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.re = tmp.re;
		out_tmp.e2.im = tmp.im;
		
		// multiply by sigma because gaussianNormalPair generates a couple
		// of real gaussian number distributed with variance 1 instead of 0.5 (see the 
		// documentation of the function set_gaussian_spinorfield_device in spinors_staggered.hpp)
		// Remark that here the variable sigma is the VARIANCE even if it is called sigma!!!
		
		out_tmp = su3vec_times_real(out_tmp, sqrt(sigma));

		put_su3vec_to_field_eo(out, id_mem, out_tmp);
	}

	prng_storeState(rngStates, &rnd);
}
