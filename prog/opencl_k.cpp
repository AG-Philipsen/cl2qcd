#include "opencl_k.h"

using namespace std;

void Opencl_k::fill_collect_options(stringstream* collect_options)
{

	Opencl_heatbath::fill_collect_options(collect_options);
	return;
}

void Opencl_k::fill_buffers()
{
	Opencl_heatbath::fill_buffers();


	cl_int clerr = CL_SUCCESS;

	cl_uint numgrps;
	size_t local_work_size;
	size_t global_work_size;
	get_work_sizes(&local_work_size, &global_work_size, &numgrps, Opencl::get_device_type() );

	cout << "Create buffer for transport coefficient kappa_karsch..." << endl;
	clmem_kappa_karsch = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);

	cout << "Create buffer for transport coefficient kappa_clover..." << endl;
	clmem_kappa_clover = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);

	return;
}

void Opencl_k::fill_kernels()
{
	Opencl_heatbath::fill_kernels();

	cout << "Create TK kappa kernels..." << endl;
	kappa_karsch_gpu = createKernel("kappa_karsch_gpu") << basic_opencl_code << "opencl_tk_kappa.cl";

	cout << "Create TK clover kernels..." << endl;
	kappa_clover_gpu = createKernel("kappa_clover_gpu") << basic_opencl_code << "opencl_tk_kappa.cl";

	return;
}



void Opencl_k::run_kappa_karsch_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_karsch_out)
{
	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	string kernelname = "dummy";
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl::get_device_type(), kernelname);	


	//buffers for kappa
	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;

	if( clmem_kappa_karsch_buf_glob == 0 ) {
		clmem_kappa_karsch_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);
	}


	clerr = clSetKernelArg(kappa_karsch_gpu, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_karsch_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_karsch_gpu, 2, sizeof(cl_mem), &clmem_kappa_karsch_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clSetKernelArg",__FILE__,__LINE__);

	enqueueKernel(kappa_karsch_gpu, global_work_size, local_work_size);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_kappa_karsch_buf_glob, CL_FALSE, 0, sizeof(hmc_float), kappa_karsch_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clEnqueueReadBuffer",__FILE__,__LINE__);

	// wait for results to have been read back
	clFinish(queue);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);

	timer->add();
	return;

}


void Opencl_k::run_kappa_clover_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_clover_out)
{


	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	string kernelname = "dummy";
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl::get_device_type(), kernelname);	


	//buffers for kappa
	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;

	if( clmem_kappa_clover_buf_glob == 0 ) {
		clmem_kappa_clover_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);
	}

	clerr = clSetKernelArg(kappa_clover_gpu, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 2, sizeof(cl_mem), &clmem_kappa_clover_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clKernelArg",__FILE__,__LINE__);

	enqueueKernel(kappa_clover_gpu, global_work_size, local_work_size);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_kappa_clover_buf_glob, CL_FALSE, 0, sizeof(hmc_float), kappa_clover_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clEnqueueReadBuffer",__FILE__,__LINE__);

	// wait for results to have been read back
	clFinish(queue);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);


	timer->add();
	return;

}

void Opencl_k::finalize_k()
{

  cl_int clerr = CL_SUCCESS;

  //	if(get_init_status() == 1) {
	  /*
		clerr = clFlush(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFlush",__FILE__,__LINE__);
		clerr = clFinish(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
	  */

	  /* //this should be part of opencl_heatbath
		clerr = clReleaseKernel(heatbath_even);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(heatbath_odd);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

		clerr = clReleaseKernel(overrelax_even);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(overrelax_odd);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

		clerr = clReleaseKernel(plaquette);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(polyakov);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(plaquette_reduction);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(polyakov_reduction);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
	  */

		clerr = clReleaseKernel(kappa_karsch_gpu);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);
		clerr = clReleaseKernel(kappa_clover_gpu);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseKernel",__FILE__,__LINE__);

		/*
		clerr = clReleaseMemObject(clmem_gaugefield);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		clerr = clReleaseMemObject(clmem_rndarray);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		*/

		/*
		clerr = clReleaseMemObject(clmem_plaq);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		clerr = clReleaseMemObject(clmem_tplaq);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		clerr = clReleaseMemObject(clmem_splaq);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		clerr = clReleaseMemObject(clmem_polyakov);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		*/

		clerr = clReleaseMemObject(clmem_kappa_karsch);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		clerr = clReleaseMemObject(clmem_kappa_clover);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);

		/*
		if(clmem_plaq_buf_glob) {
		  clerr = clReleaseMemObject(clmem_plaq_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}
		if(clmem_tplaq_buf_glob){
		  clerr = clReleaseMemObject(clmem_tplaq_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}
		if(clmem_splaq_buf_glob){
		  clerr = clReleaseMemObject(clmem_splaq_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}
		if(clmem_polyakov_buf_glob){
		  clerr = clReleaseMemObject(clmem_polyakov_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}
		*/

		if(clmem_kappa_karsch_buf_glob){
		  clerr = clReleaseMemObject(clmem_kappa_karsch_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}
		if(clmem_kappa_clover_buf_glob){ 
		  clerr = clReleaseMemObject(clmem_kappa_clover_buf_glob);
		  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);
		}

		/*
		clerr = clReleaseCommandQueue(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseCommandQueue",__FILE__,__LINE__);
		clerr = clReleaseContext(context);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseContext",__FILE__,__LINE__);
		*/

		//		set_init_false();
		//	}
	return;
}
