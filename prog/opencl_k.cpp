#include "opencl_k.h"

using namespace std;

hmc_error Opencl_k::fill_kernels_file ()
{
	//give a list of all kernel-files
	hmc_error err = Opencl::fill_kernels_file();
	if(err != HMC_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	cl_kernels_file.push_back("opencl_tk_kappa.cl");
	return HMC_SUCCESS;
}

hmc_error Opencl_k::fill_collect_options(stringstream* collect_options)
{

	hmc_error err = Opencl::fill_collect_options(collect_options);
	if(err != HMC_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	return HMC_SUCCESS;
}

hmc_error Opencl_k::fill_buffers()
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

hmc_error Opencl_k::fill_kernels()
{
	hmc_error err = Opencl::fill_kernels();

	if(err != HMC_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

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



hmc_error Opencl_k::run_kappa_karsch_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_karsch_out)
{
	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sens on CPU

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
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


hmc_error Opencl_k::run_kappa_clover_gpu(const hmc_float beta, usetimer * timer, hmc_float * kappa_clover_out)
{


	//variables
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sens on CPU

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
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

hmc_error Opencl_k::finalize()
{
	if(get_init_status() == 1) {
		if(clFlush(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clFinish(queue) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(heatbath_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(heatbath_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(overrelax_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(overrelax_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(plaquette) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(plaquette_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(polyakov_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(kappa_karsch_gpu) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(kappa_clover_gpu) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseProgram(clprogram) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseMemObject(clmem_gaugefield) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_rndarray) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseMemObject(clmem_plaq) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_tplaq) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_splaq) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseMemObject(clmem_kappa_karsch) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_kappa_clover) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clmem_plaq_buf_glob) if(clReleaseMemObject(clmem_plaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clmem_tplaq_buf_glob) if(clReleaseMemObject(clmem_tplaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clmem_splaq_buf_glob) if(clReleaseMemObject(clmem_splaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clmem_polyakov_buf_glob) if(clReleaseMemObject(clmem_polyakov_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clmem_kappa_karsch_buf_glob) if(clReleaseMemObject(clmem_kappa_karsch_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clmem_kappa_clover_buf_glob) if(clReleaseMemObject(clmem_kappa_clover_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseCommandQueue(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseContext(context) != CL_SUCCESS) exit(HMC_OCLERROR);

		set_init_false();
	}
	return HMC_SUCCESS;
}
