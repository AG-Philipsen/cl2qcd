#include "opencl_hmc.h"

hmc_error Opencl_hmc::fill_kernels_file ()
{
	//give a list of all kernel-files
	Opencl_fermions::fill_kernels_file();

	cl_kernels_file.push_back("types_hmc.h");
	cl_kernels_file.push_back("operations_gaugemomentum.cl");
	cl_kernels_file.push_back("operations_force.cl");
// 	cl_kernels_file.push_back("integrator.cl");
// 	cl_kernels_file.push_back("molecular_dynamics.cl");

	return HMC_SUCCESS;  
}

hmc_error Opencl_hmc::fill_collect_options(stringstream* collect_options)
{

	Opencl_fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << GAUGEMOMENTASIZE2;
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fill_buffers()
{
	Opencl_fermions::fill_buffers();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fill_kernels()
{
	Opencl_fermions::fill_kernels();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters)
{
	hmc_error err = Opencl_fermions::init(wanted_device_type, timer, parameters);
	err |= init_fermion_variables(parameters, timer);
	return err;
}

hmc_error Opencl_hmc::init_hmc_variables(inputparameters* parameters, usetimer * timer)
{
	(*timer).reset();

	/** @todo insert variables needed */


	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::finalize_hmc(){
	

  return HMC_SUCCESS;
}

