#include "opencl_module_heatbath.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Heatbath::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Ran::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta();
	if(get_parameters()->get_use_aniso() == true) {
		*collect_options << " -D_ANISO_";
		*collect_options <<  " -DXI_0=" << get_parameters()->get_xi_0();
	}

	return;
}


void Opencl_Module_Heatbath::fill_buffers()
{

	Opencl_Module_Ran::fill_buffers();

	return;
}

void Opencl_Module_Heatbath::fill_kernels()
{
	Opencl_Module_Ran::fill_kernels();

	logger.debug() << "Create heatbath kernels...";
	heatbath_even = createKernel("heatbath_even") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "heatbath_even.cl";
	heatbath_odd = createKernel("heatbath_odd") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "heatbath_odd.cl";
	heatbath_even_hack = createKernel("heatbath_even_hack") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "heatbath_even.cl";
	heatbath_odd_hack = createKernel("heatbath_odd_hack") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "heatbath_odd.cl";

	logger.debug() << "Create overrelax kernels...";
	overrelax_even = createKernel("overrelax_even") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "overrelax_even.cl";
	overrelax_odd = createKernel("overrelax_odd") << basic_opencl_code << prng_code << "operations_heatbath.cl" << "overrelax_odd.cl";

	return;
}

void Opencl_Module_Heatbath::clear_kernels()
{
	Opencl_Module_Ran::clear_kernels();

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(heatbath_even);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(heatbath_odd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(heatbath_even_hack);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(heatbath_odd_hack);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(overrelax_even);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(overrelax_odd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	return;
}

void Opencl_Module_Heatbath::clear_buffers()
{
	Opencl_Module_Ran::clear_buffers();


	return;
}

void Opencl_Module_Heatbath::run_heatbath()
{
	cl_int clerr = CL_SUCCESS;

	cl_mem tmp = create_rw_buffer(getGaugefieldBufferSize());
	cl_mem src = get_gaugefield();

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(heatbath_even, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_even, 1, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_even, 3, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_even_hack, 0, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_even_hack, 1, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_even, global_work_size, ls);

		clerr = clSetKernelArg(heatbath_even_hack, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_even_hack, global_work_size, ls);
	}

	this->get_work_sizes(heatbath_odd, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_odd, 1, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_odd, 3, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_odd_hack, 0, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_odd_hack, 1, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_odd, global_work_size, ls);

		clerr = clSetKernelArg(heatbath_odd_hack, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_odd_hack, global_work_size, ls);
	}

	// wait for kernel to finish to avoid hangups on AMD
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clReleaseMemObject(tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
}

void Opencl_Module_Heatbath::run_overrelax()
{
	cl_int clerr = CL_SUCCESS;

	cl_mem tmp = create_rw_buffer(getGaugefieldBufferSize());
	cl_mem src = get_gaugefield();

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(overrelax_even, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_even, 1, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_even, 3, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_even_hack, 0, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_even_hack, 1, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(overrelax_even, global_work_size, ls);

		clerr = clSetKernelArg(heatbath_even_hack, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_even_hack, global_work_size, ls);
	}

	this->get_work_sizes(overrelax_odd, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_odd, 1, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_odd, 3, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_odd_hack, 0, sizeof(cl_mem), &src);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_odd_hack, 1, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(overrelax_odd, global_work_size, ls);

		clerr = clSetKernelArg(heatbath_odd_hack, 2, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_odd_hack, global_work_size, ls);
	}

	// wait for kernel to finish to avoid hangups on AMD
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clReleaseMemObject(tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
}

void Opencl_Module_Heatbath::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Ran::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	//Query kernel name
	string kernelname = get_kernel_name(kernel);

	//Query specific sizes for kernels if needed
	//all of the following kernels are called with EnqueueKernel(gs), ls, num_groups are not needed!
	if (kernelname.compare("heatbath_even") == 0 || kernelname.compare("heatbath_odd") == 0 || kernelname.compare("overrelax_even") == 0 || kernelname.compare("overrelax_odd") == 0) {
		if( get_device_type() == CL_DEVICE_TYPE_GPU ) {
			*gs = min(parameters->get_volspace() * parameters->get_nt() / 2, this->Opencl_Module_Ran::get_num_rndstates());
		} else {
			*gs = min(get_max_compute_units(), this->Opencl_Module_Ran::get_num_rndstates());
		}
		*ls = Opencl_Module::get_numthreads();
		*num_groups = *gs / *ls;
	}
	return;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Heatbath::get_timer(const char * in)
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Ran::get_timer(in);
	if(noop != NULL) return noop;
	if (strcmp(in, "heatbath_even") == 0) {
		return &this->timer_heatbath_even;
	}
	if (strcmp(in, "heatbath_odd") == 0) {
		return &this->timer_heatbath_odd;
	}
	if (strcmp(in, "overrelax_even") == 0) {
		return &this->timer_overrelax_even;
	}
	if (strcmp(in, "overrelax_odd") == 0) {
		return &this->timer_overrelax_odd;
	}
	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}
#endif

int Opencl_Module_Heatbath::get_read_write_size(const char * in)
{
	int result = Opencl_Module_Ran::get_read_write_size(in);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	//factor for complex numbers
	int C = 2;
	const size_t VOL4D = parameters->get_vol4d();
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = get_parameters()->get_spinorfieldsize();
	//this is the same as in the function above
	if ( (strcmp(in, "heatbath_even") == 0 ) || (strcmp(in, "heatbath_odd") == 0) || (strcmp(in, "overrelax_even") == 0) || (strcmp(in, "overrelax_odd") == 0)) {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 su3-matrix
		return VOL4D / 2 * C * D * R * (6 * (NDIM - 1) + 1 + 1 );
	}
	return 0;
}

int Opencl_Module_Heatbath::get_flop_size(const char * in)
{
	int result = Opencl_Module_Ran::get_flop_size(in);
	if (result != 0) return result;
	const size_t VOL4D = parameters->get_vol4d();
	int S;
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = get_parameters()->get_spinorfieldsize();
	//this is the same as in the function above
	///@NOTE: I do not distinguish between su3 and 3x3 matrices. This is a difference if one use e.g. REC12, but here one wants to have the "netto" flops for comparability.
	if ( (strcmp(in, "heatbath_even") == 0 ) || (strcmp(in, "heatbath_odd") == 0) ) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 80 flops for the su2 update)
		return VOL4D * (4 * (NDIM - 1) * get_parameters()->get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + NC * (2 * get_parameters()->get_flop_su3_su3() + 80));
	}
	if ( (strcmp(in, "overrelax_even") == 0) || (strcmp(in, "overrelax_odd") == 0)) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 58 flops for the su2 update)
		return VOL4D * (4 * (NDIM - 1) * get_parameters()->get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + NC * (2 * get_parameters()->get_flop_su3_su3() + 58));
	}
	return 0;
}

#ifdef _PROFILING_

void Opencl_Module_Heatbath::print_profiling(std::string filename, int number)
{
	Opencl_Module_Ran::print_profiling(filename, number);
	const char * kernelName;
	kernelName = "heatbath_even";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "heatbath_odd";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "overrelax_even";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
	kernelName = "overrelax_odd";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName), this->get_flop_size(kernelName) );
}

#endif
