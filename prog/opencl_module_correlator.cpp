#include "opencl_module_correlator.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module_Correlator::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Spinors::fill_collect_options(collect_options);
	return;
}


void Opencl_Module_Correlator::fill_buffers()
{

	Opencl_Module_Spinors::fill_buffers();

	return;
}

void Opencl_Module_Correlator::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();
	return;
}

void Opencl_Module_Correlator::clear_kernels()
{
	Opencl_Module_Spinors::clear_kernels();
	return;
}

void Opencl_Module_Correlator::clear_buffers()
{
	Opencl_Module_Spinors::clear_buffers();

	return;
}

void Opencl_Module_Correlator::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
  Opencl_Module_Spinors::get_work_sizes(kernel, dev_type, ls, gs, num_groups);
  
	return;
}
