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
	return HMC_SUCCESS;  
}

hmc_error Opencl_k::fill_kernels(){
  	hmc_error err = Opencl::fill_kernels();
	return HMC_SUCCESS;  
}
