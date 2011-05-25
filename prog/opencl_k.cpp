#include "opencl_k.h"

using namespace std;

hmc_error Opencl_k::fill_kernels_file (){
	//give a list of all kernel-files
	hmc_error err =	Opencl::fill_kernels_file();
	cl_kernels_file.push_back("opencl_tk_kappa.cl");
	return HMC_SUCCESS;  
}

hmc_error Opencl_k::fill_collect_options(stringstream* collect_options){

	hmc_error err = Opencl::fill_collect_options(collect_options);
	return HMC_SUCCESS;  
}

hmc_error Opencl_k::fill_buffers(){
	hmc_error err = Opencl::fill_buffers();
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

hmc_error Opencl_k::fill_kernels(){
  	hmc_error err = Opencl::fill_kernels();
	cl_int clerr = CL_SUCCESS;
	
	cout << "Create TK kappa kernels..." << endl;
	kappa_karsch_gpu = clCreateKernel(clprogram, "kappa_karsch_gpu", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	
	cout << "Create TK clover kernels..." << endl;
	kappa_clover_gpu = clCreateKernel(clprogram, "kappa_clover_gpu", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	
	
	return HMC_SUCCESS;  
}



hmc_error Opencl_k::run_kappa_karsch_gpu(const hmc_float beta, usetimer * const timer, hmc_float * kappa_karsch_out)
{
	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();
	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	
#ifdef _USE_GPU_
	size_t global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
#else
	size_t global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);
#endif
	
	global_work_size = local_work_size * num_groups;
	
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


hmc_error Opencl_k::run_kappa_clover_gpu(const hmc_float beta, usetimer * const timer, hmc_float * kappa_clover_out)
{
	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();
	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	
#ifdef _USE_GPU_
	size_t global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
#else
	size_t global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);
#endif
	
	global_work_size = local_work_size * num_groups;
	
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


