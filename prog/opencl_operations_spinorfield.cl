//opencl_operations_spinorfield

//eoprec operations
void convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in){
  int spacepos, timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
				get_even_site(n, &spacepos, &timepos);
				even[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
				get_odd_site(n, &spacepos, &timepos);
				odd[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
      }
    }
  }
  return;
}

void convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out){
  int spacepos, timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
				get_even_site(n, &spacepos, &timepos);
				out[spinor_field_element(alpha,color,spacepos,timepos)] = even[eoprec_spinor_field_element(alpha,color,n)];
				get_odd_site(n, &spacepos, &timepos);
				out[spinor_field_element(alpha,color,spacepos,timepos)] = odd[eoprec_spinor_field_element(alpha,color,n)];
      }
    }
  }
  return;
}

void get_spinor_from_eoprec_field(__global hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[eoprec_spinor_field_element(alpha,color,n_eoprec)];
    }
  }
  return;
}

void put_spinor_to_eoprec_field(hmc_spinor* in,__global hmc_eoprec_spinor_field* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[eoprec_spinor_field_element(alpha,color,n_eoprec)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}

void get_spinor_from_field(__global hmc_spinor_field* in, hmc_spinor* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[spinor_field_element(alpha,color,n,t)];
    }
  }
  return;
}

void put_spinor_to_field(hmc_spinor* in, __global hmc_spinor_field* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_field_element(alpha,color,n,t)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}

// -alpha*x + y
//CP: defined with a minus!!!
__kernel void saxpy(__global hmc_spinor_field* x, __global hmc_spinor_field* y, __global hmc_complex * alpha, __global hmc_spinor_field* out){
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	//!!CP: Old and wrong, but perhaps one can correct this...
// 	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
		//!!CP: complexmult cannot handle __global
		hmc_complex tmp2 = x[id_tmp];
		hmc_complex tmp1 = complexmult(&alpha_tmp,&tmp2);
		((out)[id_tmp]).re = -(tmp1).re + y[id_tmp].re;
		((out)[id_tmp]).im = -(tmp1).im + y[id_tmp].im;
	}

	return;
}

//alpha*x + beta*y + z
__kernel void saxsbypz(__global hmc_spinor_field* x, __global hmc_spinor_field* y, __global hmc_spinor_field* z, __global hmc_complex * alpha, __global hmc_complex * beta, __global hmc_spinor_field* out){
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	hmc_complex beta_tmp = (*beta);
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
		//!!CP: complexmult cannot handle __global
		hmc_complex tmp2 = x[id_tmp];
		hmc_complex tmp1 = complexmult(&alpha_tmp,&tmp2);

		hmc_complex tmp3 = y[id_tmp];
		hmc_complex tmp4 = complexmult(&beta_tmp,&tmp3);		
		
    ((out)[id_tmp]).re = (tmp1).re + (tmp4).re + z[id_tmp].re;
    ((out)[id_tmp]).im = (tmp1).im + (tmp4).im + z[id_tmp].im;
	}

	return;
}

// complex (!!!) scalarproduct, return in result 
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void scalar_product( __global hmc_spinor_field *x, __global hmc_spinor_field *y, __global hmc_complex* result, __local hmc_complex* result_local ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);
	
	hmc_complex sum;
	sum.re = 0.;
	sum.im = 0.;

	//!! CP: perhaps here one can first copy a whole spinor and then do the dot-prod
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
    sum.re += x[id_tmp].re*y[id_tmp].re + x[id_tmp].im*y[id_tmp].im;
    sum.im += x[id_tmp].re*y[id_tmp].im - x[id_tmp].im*y[id_tmp].re;
	}
	
	
	if(local_size ==1){
		result[ group_id ].re = sum.re;
		result[ group_id ].im = sum.im;
	}
	else{
	// sync threads
	//!!CP: is this right??
	barrier(CLK_LOCAL_MEM_FENCE);
	//reduction
	if (idx<64){
		(result_local[idx]).re=sum.re;
		(result_local[idx]).im=sum.im;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=64){
		result_local[idx%64].re +=result_local[idx].re;
		result_local[idx%64].im +=result_local[idx].im;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=32){
		result_local[idx-32].re +=result_local[idx].re;
		result_local[idx-32].im +=result_local[idx].im;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=16){
		result_local[idx-16].re +=result_local[idx].re;
		result_local[idx-16].im +=result_local[idx].im;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=8){
		result_local[idx-8].re +=result_local[idx].re;
		result_local[idx-8].im +=result_local[idx].im;
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	//thread 0 sums up the result_local and stores it in array result
	if (idx==0)
		for (int i = 1; i<8; i++) {
			result_local[idx].re += result_local[i].re;
			result_local[idx].im += result_local[i].im;
// 			complexaccumulate(result[get_group_id(0)], result_local[i]);
		}
		result[ group_id ].re = result_local[idx].re;
		result[ group_id ].im = result_local[idx].im;	
	}
	return;
}

__kernel void scalar_product_reduction(__global hmc_complex* result_tmp, __global hmc_complex* result){
	//!!CP: complex_acc cannot handle __global
	hmc_complex tmp1;
	hmc_complex tmp2;
	tmp2 = result_tmp[0];
	int id = get_global_id(0);
	//!!CP: only one thread is needed
	if(id == 0){
		for (int i = 1; i<get_num_groups(0); i++){
			tmp1 = result_tmp[i];
			complexaccumulate(&tmp2, &tmp1);
		}
	//!!CP: pointer????
		(*result) = tmp2;
	}
	return;
}

// complex (!!!) squarenorm, return in result 
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults 
__kernel void global_squarenorm( __global hmc_spinor_field *x, __global hmc_float* result, __local hmc_float* result_local ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);
	
	hmc_float sum;
	sum = 0.;

	//!! CP: perhaps here one can first copy a whole spinor and then do the dot-prod
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
    sum += x[id_tmp].re*x[id_tmp].re + x[id_tmp].im*x[id_tmp].im;
	}
	
	if(local_size ==1){
		result[ group_id ] = sum;
	}
	else{
	// sync threads
	//!!CP: is this right, see also above??
	barrier(CLK_LOCAL_MEM_FENCE);
	//reduction
	if (idx<64)
		(result_local[idx])=sum;
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=64)
		result_local[idx%64]+=result_local[idx];
//!! CP: this must be result_local[id%64]+=sum instead of [id-64]; ??
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=32)
		result_local[idx-32]+=result_local[idx];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=16)
		result_local[idx-16]+=result_local[idx];
	barrier(CLK_LOCAL_MEM_FENCE);
	if (idx>=8)
		result_local[idx-8]+=result_local[idx];
	barrier(CLK_LOCAL_MEM_FENCE);
	//thread 0 sums up the result_local and stores it in array result
	if (idx==0)
		result[ group_id ] = 	result_local[0] + result_local[1] + result_local[2] + result_local[3] + 
				                result_local[4] + result_local[5] + result_local[6] + result_local[7];
	}
	return;
}

__kernel void global_squarenorm_reduction(__global hmc_float* result_tmp, __global hmc_float* result){
	int id = get_global_id(0);
	//!!CP: only one thread is needed
	if(id == 0){
		for (int i = 1; i<get_num_groups(0); i++){
			result_tmp[0] += result_tmp[i];
		}
		//!!CP: pointer????
		(*result) = result_tmp[0];
	}

	return;
}

__kernel void set_zero_spinorfield( __global hmc_spinor_field *x ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
		x[id_tmp].re = 0.;
		x[id_tmp].im = 0.;
	}
	return;
}

__kernel void convert_to_kappa_format( __global hmc_spinor_field *in, __global hmc_float * kappa ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_float tmp = *kappa;
	
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
		in[id_tmp].re *= sqrt(2.*tmp);
		in[id_tmp].im *= sqrt(2.*tmp);
	}
	return;
}

__kernel void convert_from_kappa_format( __global hmc_spinor_field *in, __global hmc_spinor_field* out, __global hmc_float * kappa ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_float tmp = *kappa;
	
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size){
		out[id_tmp].re = (in[id_tmp].re)/sqrt(2.*tmp);
		out[id_tmp].im = (in[id_tmp].im)/sqrt(2.*tmp);
	}
	return;
}

__kernel void convert_to_kappa_format_eoprec( __global hmc_spinor_field *in, __global hmc_float * kappa ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_float tmp = *kappa;
	
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size){
		in[id_tmp].re *= sqrt(2.*tmp);
		in[id_tmp].im *= sqrt(2.*tmp);
	}
	return;
}

__kernel void convert_from_kappa_format_eoprec( __global hmc_spinor_field *in, __global hmc_spinor_field* out, __global hmc_float * kappa ){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	
	hmc_float tmp = *kappa;
	
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size){
		out[id_tmp].re = (in[id_tmp].re)/sqrt(2.*tmp);
		out[id_tmp].im = (in[id_tmp].im)/sqrt(2.*tmp);
	}
	return;
}

//!!CP: these two need to be kernels, if they are needed at all...
void create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, __global hmc_ocl_gaugefield * gaugefield){

	/*
	set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

	*/
	
  return;
}

void create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta,hmc_float chem_pot_re, hmc_float chem_pot_im, __global hmc_ocl_gaugefield* gaugefield){
/*  
  see host code
*/
  return;
}
