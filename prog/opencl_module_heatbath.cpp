#include "opencl_module_heatbath.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Heatbath::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Ran::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta();

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
	heatbath_even = createKernel("heatbath_even") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
	heatbath_odd = createKernel("heatbath_odd") << basic_opencl_code << "random.cl" << "update_heatbath.cl";

	logger.debug() << "Create overrelax kernels...";
	overrelax_even = createKernel("overrelax_even") << basic_opencl_code << "random.cl" << "overrelax.cl";
	overrelax_odd = createKernel("overrelax_odd") << basic_opencl_code << "random.cl" << "overrelax.cl";

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

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(heatbath_even, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), get_gaugefield());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_even, 2, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel(heatbath_even, global_work_size);
	}

	this->get_work_sizes(heatbath_odd, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), get_gaugefield());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(heatbath_odd, 2, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		enqueueKernel(heatbath_odd, global_work_size);
	}

	// do not wait for kernel to finish...
	//  clerr = clFinish(get_queue());
	//  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
	return;

}

void Opencl_Module_Heatbath::run_overrelax()
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(overrelax_even, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), get_gaugefield());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(overrelax_even, 2, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel(overrelax_even, global_work_size);
	}

	this->get_work_sizes(overrelax_odd, this->get_device_type(), &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), get_gaugefield());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(overrelax_odd, 2, sizeof(cl_mem), get_clmem_rndarray());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 1, sizeof(int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

		enqueueKernel(overrelax_odd, global_work_size);
	}

	//do not wait for kernel to finish
	//  clerr = clFinish(get_queue());
	//  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
	return;
}

void Opencl_Module_Heatbath::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Ran::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	//Query kernel name
	string kernelname = get_kernel_name(kernel);

	//Query specific sizes for kernels if needed
	//all of the following kernels are called with EnqueueKernel(gs), ls, num_groups are not needed!
	if (kernelname.compare("heatbath_even") == 0 || kernelname.compare("heatbath_odd") == 0 || kernelname.compare("overrelax_even") == 0 || kernelname.compare("overrelax_even") == 0) {
		if( get_device_type() == CL_DEVICE_TYPE_GPU )
			*gs = min(VOLSPACE * NTIME / 2, this->Opencl_Module_Ran::get_num_rndstates());
		else
			*gs = min(get_max_compute_units(), this->Opencl_Module_Ran::get_num_rndstates());
		*ls = 0;
		*num_groups = 0;
	}
	return;
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Heatbath::get_timer(char * in)
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
int Opencl_Module_Heatbath::get_read_write_size(char * in, inputparameters * parameters)
{
	Opencl_Module_Ran::get_read_write_size(in, parameters);
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
		S = get_parameters()->get_eoprec_spinorfieldsize();
	else
		S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "heatbath_even") == 0) {
		return VOL4D * D * R + 1;
	}
	if (strcmp(in, "heatbath_odd") == 0) {
		return VOL4D * D * R + 1;
	}
	if (strcmp(in, "overrelax_even") == 0) {
		return 48 * VOL4D * D * R + 1;
	}
	if (strcmp(in, "overrelax_odd") == 0) {
		return 48 * VOL4D * D * R + 1;
	}
}

void Opencl_Module_Heatbath::print_profiling(std::string filename)
{
	Opencl_Module_Ran::print_profiling(filename);
	char * kernelName;
	kernelName = "heatbath_even";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "heatbath_odd";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "overrelax_even";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "overrelax_odd";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
}

#endif
