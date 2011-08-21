#include "opencl_k_hybrid.h"

using namespace std;

hmc_error Opencl_k_hybrid::fill_collect_options(stringstream* collect_options)
{

	hmc_error err = Opencl::fill_collect_options(collect_options);
	if(err != HMC_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	return HMC_SUCCESS;
}

hmc_error Opencl_k_hybrid::fill_buffers()
{
	hmc_error err = Opencl::fill_buffers();

	if(err != HMC_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cl_int clerr = CL_SUCCESS;

	cout << "Create buffer for transport coefficient kappa_karsch..." << endl;
	clmem_kappa_karsch = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Create buffer for transport coefficient kappa_clover..." << endl;
	clmem_kappa_clover = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	return HMC_SUCCESS;
}

void Opencl_k_hybrid::fill_kernels()
{
	Opencl::fill_kernels();

	cout << "Create TK kappa kernels..." << endl;
	kappa_karsch_gpu = createKernel("kappa_karsch_gpu") << basic_opencl_code << "opencl_tk_kappa.cl";

	cout << "Create TK clover kernels..." << endl;
	kappa_clover_gpu = createKernel("kappa_clover_gpu") << basic_opencl_code << "opencl_tk_kappa.cl";
}



hmc_error Opencl_k_hybrid::run_kappa_karsch_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_karsch_out)
{
	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	char * kernelname = "dummy";
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl::get_device_type(), kernelname);	

	//buffers for kappa
	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;

	if( clmem_kappa_karsch_buf_glob == 0 ) {
		clmem_kappa_karsch_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating clmem_plaq_buf_glob failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
	}

	clerr = clSetKernelArg(kappa_karsch_gpu, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 at kappa_karsch_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(kappa_karsch_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 at kappa_karsch_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(kappa_karsch_gpu, 2, sizeof(cl_mem), &clmem_kappa_karsch_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 at kappa_karsch_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	enqueueKernel(kappa_karsch_gpu, global_work_size, local_work_size);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_kappa_karsch_buf_glob, CL_FALSE, 0, sizeof(hmc_float), kappa_karsch_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	// wait for results to have been read back
	clFinish(queue);

	timer->add();
	return HMC_SUCCESS;

}


hmc_error Opencl_k_hybrid::run_kappa_clover_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_clover_out)
{


	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	char * kernelname = "dummy";
	get_work_sizes(&local_work_size, &global_work_size, &num_groups, Opencl::get_device_type(), kernelname);	


	//buffers for kappa
	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;

	if( clmem_kappa_clover_buf_glob == 0 ) {
		clmem_kappa_clover_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating clmem_plaq_buf_glob failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
	}

	clerr = clSetKernelArg(kappa_clover_gpu, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 at kappa_clover_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(kappa_clover_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 at kappa_clover_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(kappa_clover_gpu, 2, sizeof(cl_mem), &clmem_kappa_clover_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 at kappa_clover_gpu failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	enqueueKernel(kappa_clover_gpu, global_work_size, local_work_size);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_kappa_clover_buf_glob, CL_FALSE, 0, sizeof(hmc_float), kappa_clover_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	// wait for results to have been read back
	clFinish(queue);

	timer->add();
	return HMC_SUCCESS;

}

hmc_error Opencl_k_hybrid::finalize()
{
	if(get_init_status() == 1) {
		if(clFlush(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clFinish(queue) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(kappa_karsch_gpu) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(kappa_clover_gpu) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseMemObject(clmem_kappa_karsch) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_kappa_clover) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clmem_kappa_karsch_buf_glob) if(clReleaseMemObject(clmem_kappa_karsch_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clmem_kappa_clover_buf_glob) if(clReleaseMemObject(clmem_kappa_clover_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseCommandQueue(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseContext(context) != CL_SUCCESS) exit(HMC_OCLERROR);

		set_init_false();
	}
	return HMC_SUCCESS;
}
