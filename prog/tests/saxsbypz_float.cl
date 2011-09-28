//typedefs for kernels

#define SIZE 1024

typedef float hmc_float;

typedef struct {
	hmc_float re;
	hmc_float im;
} hmc_complex;

typedef struct {
	hmc_complex e0;
	hmc_complex e1;
	hmc_complex e2;
	
} su3vec;

typedef struct {
	su3vec e0;
	su3vec e1;
	su3vec e2;
	su3vec e3;
} spinor;

//functions used in kernels

su3vec su3vec_times_complex(su3vec in, hmc_complex factor){
	su3vec tmp;
	tmp.e0.re = in.e0.re*factor.re - in.e0.im*factor.im;
	tmp.e0.im = in.e0.im*factor.re + in.e0.re*factor.im;
	tmp.e1.re = in.e1.re*factor.re - in.e1.im*factor.im;
	tmp.e1.im = in.e1.im*factor.re + in.e1.re*factor.im;
	tmp.e2.re = in.e2.re*factor.re - in.e2.im*factor.im;
	tmp.e2.im = in.e2.im*factor.re + in.e2.re*factor.im;
	return tmp;
}

spinor spinor_times_complex(spinor in, hmc_complex factor)
{
	spinor tmp;
	tmp.e0 = su3vec_times_complex(in.e0, factor);
	tmp.e1 = su3vec_times_complex(in.e1, factor);
	tmp.e2 = su3vec_times_complex(in.e2, factor);
	tmp.e3 = su3vec_times_complex(in.e3, factor);
	return tmp;
}

su3vec su3vec_acc_acc(su3vec in1, su3vec in2, su3vec in3){
	su3vec tmp;     
	tmp.e0.re = in1.e0.re + in2.e0.re + in3.e0.re;
	tmp.e0.im = in1.e0.im + in2.e0.im + in3.e0.im;
	tmp.e1.re = in1.e1.re + in2.e1.re + in3.e1.re;
	tmp.e1.im = in1.e1.im + in2.e1.im + in3.e1.im;
	tmp.e2.re = in1.e2.re + in2.e2.re + in3.e2.re;
	tmp.e2.im = in1.e2.im + in2.e2.im + in3.e2.im;
	return tmp;
}

spinor spinor_acc_acc(spinor in1, spinor in2, spinor in3)
{
	spinor tmp;
	tmp.e0 = su3vec_acc_acc(in1.e0, in2.e0, in3.e0);
	tmp.e1 = su3vec_acc_acc(in1.e1, in2.e1, in3.e1);
	tmp.e2 = su3vec_acc_acc(in1.e2, in2.e2, in3.e2);
	tmp.e3 = su3vec_acc_acc(in1.e3, in2.e3, in3.e3);
	return tmp;
}

//2 kernels to calculate
//	alpha*x + beta*y + z
//with 	alpha, beta complex numbers
//and 	x,y,z complex spinor-vector of size SIZE

__kernel void saxsbypz_1(__global spinor* x, __global spinor* y, __global spinor* z, __global hmc_complex * alpha, __global hmc_complex * beta, __global spinor* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	hmc_complex beta_tmp = (*beta);
	for(int id_tmp = id; id_tmp < SIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		spinor y_tmp = y[id_tmp];
		spinor z_tmp = z[id_tmp];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		out[id_tmp] = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
	}

	return;
}

__kernel void saxsbypz_2(__global spinor* x, __global spinor* y, __global spinor* z, __global hmc_complex * alpha, __global hmc_complex * beta, __global spinor* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	hmc_complex beta_tmp = (*beta);
	for(int id_tmp = id; id_tmp < SIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = y[id_tmp];
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = z[id_tmp];
		out[id_tmp] = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
	}

	return;
}

