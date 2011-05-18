#include "opencl_k.h"

using namespace std;

hmc_error Opencl_k::fill_kernels_file (){
	//give a list of all kernel-files
	hmc_error err =	Opencl::fill_kernels_file();
	cl_kernels_file.push_back("opencl_tk_kappa.cl");
  
	return HMC_SUCCESS;  
}
