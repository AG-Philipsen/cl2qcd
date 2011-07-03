/** @file
 * Device code for operations on the even-odd preconditioned spinor field
 */

//eoprec operations
void convert_to_eoprec(spinorfield_eoprec* even, spinorfield_eoprec* odd, spinorfield* in)
{
	int spacepos, timepos;
	for(int n=0; n<VOL4D/2; n++) {
		     get_even_site(n, &spacepos, &timepos);
		     even[n] = in[get_global_pos(spacepos,timepos)];
		     get_odd_site(n, &spacepos, &timepos);
		     odd[n] = in[get_global_pos(spacepos,timepos)];
		}
	return;
}

__kernel void convert_from_eoprec(__global spinorfield_eoprec* even, __global spinorfield_eoprec* odd, __global spinorfield* out){
  int id = get_global_id(0);
  if(id ==0){
	int spacepos, timepos;
	for(int n=0; n<EOPREC_SPINORFIELDSIZE2; n++) {
		get_even_site(n, &spacepos, &timepos);
		out[get_global_pos(spacepos,timepos)] = even[n];
		get_odd_site(n, &spacepos, &timepos);
		out[get_global_pos(spacepos,timepos)] = odd[n];
	}
   }
	return;
}

//CP: if this works it can also be done directly in the code!!
spinor get_spinor_from_eoprec_field(__global spinorfield_eoprec* in, int n_eoprec)
{
	spinor out;
	out = in[n_eoprec];
	return out;
}

void put_spinor_to_eoprec_field(spinor in,__global spinorfield_eoprec * out, int n_eoprec)
{
	out[n_eoprec] = in;
}

// -alpha*x + y
//CP: defined with a minus!!!
__kernel void saxpy_eoprec(__global spinorfield_eoprec* x, __global spinorfield_eoprec* y, __global hmc_complex * alpha, __global spinorfield_eoprec* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = y[id_tmp];
		out[id_tmp] = spinor_dim(y_tmp, x_tmp);
	}
}

//alpha*x + beta*y + z
__kernel void saxsbypz_eoprec(__global spinorfield_eoprec* x, __global spinorfield_eoprec* y, __global spinorfield_eoprec* z, __global hmc_complex * alpha, __global hmc_complex * beta, __global spinorfield_eoprec* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	hmc_complex beta_tmp = (*beta);
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = y[id_tmp];
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = z[id_tmp];
		
		out[id_tmp] = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
	}

	return;
}

// complex (!!!) scalarproduct, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void scalar_product_eoprec( __global spinorfield_eoprec *x, __global spinorfield_eoprec *y, __global hmc_complex* result, __local hmc_complex* result_local )
{
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

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		spinor y_tmp = y[id_tmp];
		hmc_complex tmp = spinor_scalarproduct(x_tmp, y_tmp);
		sum.re += tmp.re;
		sum.im += tmp.im;
	}


	if(local_size ==1) {
		result[ group_id ].re = sum.re;
		result[ group_id ].im = sum.im;
	} else {
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		(result_local[idx]).re=sum.re;
		(result_local[idx]).im=sum.im;

		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx>=64) {
			result_local[idx%64].re +=result_local[idx].re;
			result_local[idx%64].im +=result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx>=32) {
			result_local[idx-32].re +=result_local[idx].re;
			result_local[idx-32].im +=result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx>=16) {
			result_local[idx-16].re +=result_local[idx].re;
			result_local[idx-16].im +=result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx>=8) {
			result_local[idx-8].re +=result_local[idx].re;
			result_local[idx-8].im +=result_local[idx].im;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		//thread 0 sums up the result_local and stores it in array result
		if (idx==0) {
			if(local_size >= 8) {
				for (int i = 1; i<8; i++) {
					result_local[idx].re += result_local[i].re;
					result_local[idx].im += result_local[i].im;
				}
				result[ group_id ].re = result_local[idx].re;
				result[ group_id ].im = result_local[idx].im;
			} else {
				for (int i = 1; i<local_size; i++) {
					result_local[idx].re += result_local[i].re;
					result_local[idx].im += result_local[i].im;
				}
				result[ group_id ].re = result_local[idx].re;
				result[ group_id ].im = result_local[idx].im;
			}
		}
	}
	return;
}

// hmc_float squarenorm, return in result
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void global_squarenorm_eoprec( __global spinorfield_eoprec *x, __global hmc_float* result, __local hmc_float* result_local )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int idx = get_local_id(0);

	hmc_float sum;
	sum = 0.;

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		hmc_float tmp = spinor_squarenorm(x_tmp);
		sum += tmp;
	}

	if(local_size ==1) {
		result[ group_id ] = sum;
	} else {
		// sync threads
		barrier(CLK_LOCAL_MEM_FENCE);
		//reduction
		(result_local[idx])=sum;
		barrier(CLK_LOCAL_MEM_FENCE);
		if (idx>=64)
			result_local[idx%64]+=result_local[idx];
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
		if (idx==0) {
			if(local_size >= 8) {
				result[ group_id ] = 	result_local[0] + result_local[1] + result_local[2] + result_local[3] +
				                      result_local[4] + result_local[5] + result_local[6] + result_local[7];
			} else {
				for(int i = 0; i<local_size; i++)
					result[group_id] += result_local[i];
			}
		}
	}
	return;
}

__kernel void set_zero_spinorfield_eoprec( __global spinorfield_eoprec *x )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		x[id_tmp] = set_spinor_zero();
	}
	return;
}

__kernel void convert_to_kappa_format_eoprec( __global spinorfield_eoprec *in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_float tmp = sqrt(KAPPA*2.);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		in[id_tmp] = real_multiply_spinor(in[id_tmp], tmp);
	}
	return;
}

__kernel void convert_from_kappa_format_eoprec( __global spinorfield_eoprec *in, __global spinorfield_eoprec* out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_float tmp = 1./sqrt(KAPPA*2.);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		in[id_tmp] = real_multiply_spinor(in[id_tmp], tmp);
	}
	return;
}

__kernel void create_point_source_eoprec(__global spinorfield_eoprec* b, int i, int n)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_float tmp = sqrt(2.*KAPPA);
		int color = spinor_color(i);
		int spin = spinor_spin(i,color);
		int pos = n;

		switch (color){
			
			case 0:
				switch (spin){
					case 0:
						(b[pos].e0).e0.re = tmp;
						break;
					case 1:
						(b[pos].e1).e0.re = tmp;
						break;
					case 2:
						(b[pos].e2).e0.re = tmp;
						break;
					case 3:
						(b[pos].e3).e0.re = tmp;
						break;
				}
				break;
			case 1:
				switch (spin){
					case 0:
						(b[pos].e0).e1.re = tmp;
						break;
					case 1:
						(b[pos].e1).e1.re = tmp;
						break;
					case 2:
						(b[pos].e2).e1.re = tmp;
						break;
					case 3:
						(b[pos].e3).e1.re = tmp;
						break;
				}
				break;
			case 2:
				switch (spin){
					case 0:
						(b[pos].e0).e2.re = tmp;
						break;
					case 1:
						(b[pos].e1).e2.re = tmp;
						break;
					case 2:
						(b[pos].e2).e2.re = tmp;
						break;
					case 3:
						(b[pos].e3).e2.re = tmp;
						break;
				}
				break;
}
	}
	return;
}

__kernel void set_eoprec_spinorfield_cold(__global spinorfield* in){
  int id = get_global_id(0);
  if (id ==0){
     	 //NOTE: the normalization is the same as in the above case, putting |phi|^2 of the WHOLE field to 1!!
  	 hmc_float norm = 1./sqrt(12. * VOLSPACE * NTIME);
	 for(int n = 0; n<EOPREC_SPINORFIELDSIZE; n++){
	 	 in[n] = set_spinor_cold();
		 in[n] = real_multiply_spinor(in[n], norm);
	 }
  }
}

