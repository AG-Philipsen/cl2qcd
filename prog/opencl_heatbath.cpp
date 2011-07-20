#include "opencl_heatbath.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

hmc_error Opencl_heatbath::fill_collect_options(stringstream* collect_options)
{
  Opencl::fill_collect_options(collect_options);
  return HMC_SUCCESS;
}


hmc_error Opencl_heatbath::fill_buffers()
{
  Opencl::fill_buffers();

  cl_int clerr = CL_SUCCESS;
	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_rndarray), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}


	return HMC_SUCCESS;
}

hmc_error Opencl_heatbath::clear_buffers(){

  Opencl::clear_buffers();
  if(clReleaseMemObject(clmem_rndarray) != CL_SUCCESS) exit(HMC_OCLERROR);

  return HMC_SUCCESS;

}

void Opencl_heatbath::fill_kernels()
{
  Opencl::fill_kernels();

  logger.debug() << "Create heatbath kernels...";
  heatbath_even = createKernel("heatbath_even") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
  heatbath_odd = createKernel("heatbath_odd") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
  
  logger.debug() << "Create overrelax kernels...";
  overrelax_even = createKernel("overrelax_even") << basic_opencl_code << "random.cl" << "overrelax.cl";
  overrelax_odd = createKernel("overrelax_odd") << basic_opencl_code << "random.cl" << "overrelax.cl";

}

hmc_error Opencl_heatbath::clear_kernels()
{
  Opencl::clear_kernels();

		if(clReleaseKernel(heatbath_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(heatbath_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(overrelax_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(overrelax_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

	return HMC_SUCCESS;
}

hmc_error Opencl_heatbath::run_heatbath(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_even failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_even, global_work_size);
	}

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_odd failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;

}

hmc_error Opencl_heatbath::run_overrelax(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;

	timer->reset();

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg4 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_even, global_work_size);
	}

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg7 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg8 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;
}


hmc_error Opencl_heatbath::copy_rndarray_to_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Copy randomarray to device..."<<endl;
	timer->reset();

	int clerr = clEnqueueWriteBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl_heatbath::copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Get randomarray from device..."<<endl;
	timer->reset();

	int clerr = clEnqueueReadBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

