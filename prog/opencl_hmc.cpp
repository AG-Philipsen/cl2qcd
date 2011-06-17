#include "opencl_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

hmc_error Opencl_hmc::fill_kernels_file ()
{
	//give a list of all kernel-files
	Opencl_fermions::fill_kernels_file();

	cl_kernels_file.push_back("types_hmc.h");
	cl_kernels_file.push_back("operations_gaugemomentum.cl");
	cl_kernels_file.push_back("operations_force.cl");
 	cl_kernels_file.push_back("molecular_dynamics.cl");

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
	err |= init_hmc_variables(parameters, timer);
	return err;
}

hmc_error Opencl_hmc::init_hmc_variables(inputparameters* parameters, usetimer * timer)
{
	(*timer).reset();
	
	//CP: this is copied from opencl_fermions
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

	(*timer).reset();

	logger.trace()<< "init HMC variables...";
	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor)*SPINORFIELDSIZE;
	int eoprec_spinorfield_size = sizeof(spinor)*EOPREC_SPINORFIELDSIZE;
	int gaugemomentum_size = sizeof(ae)*GAUGEMOMENTASIZE2;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
	int global_buf_size = complex_size * num_groups;
	int global_buf_size_float = float_size * num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	hmc_float tmp;
	

	/** @todo insert variables needed */
	//init mem-objects

	logger.trace() << "Create buffer for force...";
	clmem_force = clCreateBuffer(context,CL_MEM_READ_WRITE,gaugemomentum_size,0,&clerr);;
	if(clerr!=CL_SUCCESS) {
		cout<<"creating clmem_force failed, aborting..."<<endl;
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for phi_inv...";
	clmem_phi_inv = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
	if(clerr!=CL_SUCCESS) {
		cout<<"creating clmem_phi_inv failed, aborting..."<<endl;
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for new gaugefield...";
	clmem_new_u = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(s_gaugefield), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for gaugemomentum p...";
	clmem_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for gaugemomentum new_p...";
	clmem_new_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for initial spinorfield energy...";
	clmem_energy_init = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "Create buffer for deltaH...";
	clmem_deltah = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	//init kernels
	logger.debug() << "Create kernel set_zero_gaugemomentum...";
	set_zero_gaugemomentum = clCreateKernel(clprogram, "set_zero_gaugemomentum", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel generate_gaussian_spinorfield...";
	generate_gaussian_spinorfield = clCreateKernel(clprogram, "generate_gaussian_spinorfield", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel generate_gaussian_gaugemomenta...";
	generate_gaussian_gaugemomenta = clCreateKernel(clprogram, "generate_gaussian_gaugemomenta", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel md_update_gaugefield...";
	md_update_gaugefield = clCreateKernel(clprogram, "md_update_gaugefield", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel md_update_gaugemomenta...";
	md_update_gaugemomenta = clCreateKernel(clprogram, "md_update_gaugemomenta", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel gauge_force...";
	gauge_force = clCreateKernel(clprogram, "gauge_force", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	logger.debug() << "Create kernel fermion_force...";
	fermion_force = clCreateKernel(clprogram, "fermion_force", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	s_gauge = clCreateKernel(clprogram, "s_gauge", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	/** @todo this can most likely be deleted */
	s_fermion = clCreateKernel(clprogram, "s_fermion", &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}


	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::finalize_hmc(){

	logger.debug() << "release HMC-variables.." ;
	if(clReleaseMemObject(clmem_energy_init)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_deltah)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_u)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_phi_inv)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_force)!=CL_SUCCESS) return HMC_RELEASEVARIABLEERR;

	logger.debug() << "release HMC-kernels.." ;
	if(clReleaseKernel(generate_gaussian_spinorfield)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(s_gauge)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(s_fermion)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(generate_gaussian_gaugemomenta)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugefield)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugemomenta)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(gauge_force)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(fermion_force)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(set_zero_gaugemomentum)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;

	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

hmc_error Opencl_hmc::generate_gaussian_gaugemomenta_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();
	int clerr;
	//this is always applied to clmem_new_p
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta,0,sizeof(cl_mem),&clmem_new_p);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta,1,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,generate_gaussian_gaugemomenta,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue generate_gaussian_gaugemomenta kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::generate_gaussian_spinorfield_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();
	int clerr;
	//this is always applied to clmem_phi_inv
	clerr = clSetKernelArg(generate_gaussian_spinorfield,0,sizeof(cl_mem),&clmem_phi_inv);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(generate_gaussian_spinorfield,1,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,generate_gaussian_spinorfield,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue generate_gaussian_spinorfield kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::md_update_spinorfield_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	//suppose the initial gaussian field is saved in phi_inv. then the "phi" from the algorithm is clmem_inout
	int err =  Opencl_fermions::Qplus_device(clmem_phi_inv, get_clmem_inout() , local_work_size, global_work_size,  timer);

	(*timer).add();

	if(err!=HMC_SUCCESS)
		logger.fatal() << "error occured in md_update_spinorfield_device.. ";
	return HMC_SUCCESS;
}


//CP: steps2 is not used at the moment....
hmc_error Opencl_hmc::leapfrog_device(hmc_float tau, int steps1, int steps2, const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();
	//it is assumed that gaugefield and gaugemomentum have been set to the old ones already
	hmc_float stepsize = (tau) /((hmc_float) steps1);
	int k;
	hmc_float stepsize_half = 0.5*stepsize;

	//initial step
	logger.trace()<< "\t\tinitial step:";
	//here, phi is inverted using the orig. gaugefield
	/** @todo replace noop here */
	usetimer noop;
	force_device(ls, gs, &noop);
	md_update_gaugemomentum_device(stepsize_half, ls, gs, &noop);
	//intermediate steps
	if(steps1 > 1) logger.trace()<< "\t\tperform " << steps1-1 << " intermediate steps " ;
	for(k = 1; k<steps1; k++){
		md_update_gaugefield_device(stepsize, ls, gs, &noop);
		force_device(ls, gs, &noop);
		md_update_gaugemomentum_device(stepsize, ls, gs, &noop);
	}
	
	//final step
	logger.trace()<< "\t\tfinal step" ;
	md_update_gaugefield_device(stepsize, ls, gs, &noop);
	force_device(ls, gs, &noop);
	md_update_gaugemomentum_device(stepsize_half, ls, gs, &noop); 
	
	logger.trace() << "\t\tfinished leapfrog";
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::force_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();
	/** @todo replace noop */
	usetimer noop;
	cout << "\t\tstart calculating the force..." << endl;
	//CP: make sure that the output field is set to zero
	set_zero_clmem_force_device(ls,gs, &noop);
	//add contributions
	cout << "\t\tcalc gauge_force..." << endl;
	gauge_force_device(ls, gs, &noop);
	cout << "\t\tinvert fermion field..." << endl;
	//CP: to begin with, consider only the cg-solver
	//source is at 0
	int k = 0;
	int use_cg = TRUE;
	//CP: at the moment, use_eo = 0 so that even-odd is not used!!!!!
	
	//debugging
	int err = 0;
	int use_eo = 0;
	/** @todo check the use of the sources again, compare to tmlqcd!!! */
	if(use_cg){
		if(!use_eo){
			//the inversion calculates Y = (QplusQminus)^-1 phi = phi_inv
			cout << "\t\t\tstart solver" << endl;
			Opencl_fermions::create_point_source_device(k, 0,0,ls, gs, &noop);
			Opencl_fermions::solver_device(&noop,&noop,&noop,&noop,&noop,&noop,&noop,&noop, ls, gs, get_parameters()->get_cgmax());
	
			if (err != HMC_SUCCESS) cout << "\t\tsolver did not solve!!" << endl;
			else cout << "\t\tsolver solved!" << endl;
		}
		else{
			hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
			hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
			
			Opencl_fermions::create_point_source_eoprec_device(k,0,0, ls, gs, &noop, &noop, &noop);
			Opencl_fermions::solver_eoprec_device(&noop,&noop,&noop,&noop,&noop,&noop,&noop,&noop, ls, gs, get_parameters()->get_cgmax());
		}
		cout << "\t\t\tcalc X" << endl;
		//X = Qminus Y = Qminus phi_inv 
		Qminus_device(clmem_phi_inv,get_clmem_inout(), ls, gs, &noop);
	}
	else{
		//here, one has first to invert (with BiCGStab) Qplus phi = X and then invert Qminus X => Qminus^-1 Qplus^-1 phi = (QplusQminus)^-1 phi = Y = phi_inv
	}
	cout << "\t\tcalc fermion_force..." << endl;
	fermion_force_device(ls, gs, &noop);
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::hamiltonian_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	//Without fermions here!!! H = S_gauge + S_gaugemomenta
// hmc_float hamiltonian(hmc_gaugefield * field, hmc_float beta, hmc_algebraelement2 * p){
// 	hmc_float result;
// 	(result) = 0.;
// 	(result) += s_gauge(field, beta);
// 	//s_gm = 1/2*squarenorm(Pl)
// 	hmc_float s_gm;
// 	gaugemomenta_squarenorm(p, &s_gm);
// 	result += 0.5*s_gm;
// 	
// 	return result;
// }
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::calc_spinorfield_init_energy_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugemomentum_device(hmc_float eps, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugefield_device(hmc_float eps, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::set_zero_clmem_force_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();

	int clerr;
	//this is always applied to clmem_force
	clerr = clSetKernelArg(set_zero_gaugemomentum,0,sizeof(cl_mem),&clmem_force);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,set_zero_gaugemomentum,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue set_zero_gaugemomentum kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::gauge_force_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();

		int clerr;
	clerr = clSetKernelArg(gauge_force,0,sizeof(cl_mem),&clmem_new_u);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(gauge_force,1,sizeof(cl_mem),&clmem_force);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,gauge_force,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue gauge_force kernel failed, aborting..." << clerr <<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);	
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fermion_force_device(const size_t ls, const size_t gs, usetimer * timer){
	(*timer).reset();
	
	//fermion_force(field, Y, X, out);
	cl_mem tmp = get_clmem_inout();
	int clerr;
	clerr = clSetKernelArg(fermion_force,0,sizeof(cl_mem),&clmem_new_u);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(fermion_force,1,sizeof(cl_mem),&clmem_phi_inv);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(fermion_force,2,sizeof(cl_mem),&tmp);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(fermion_force,3,sizeof(cl_mem),&clmem_force);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,fermion_force,1,0,&gs,&ls,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue fermion_force kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }		
	clFinish(queue);	
	
	(*timer).add();
	return HMC_SUCCESS;

}

////////////////////////////////////////////////////
//Methods to copy new and old fields... these can be optimized!!
hmc_error Opencl_hmc::copy_gaugefield_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	int clerr = clEnqueueCopyBuffer(queue, clmem_new_u, clmem_gaugefield, 0, 0, sizeof(s_gaugefield), NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy old to new gaugefield failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	int gaugemomentum_size = sizeof(ae)*GAUGEMOMENTASIZE2;
	int clerr = clEnqueueCopyBuffer(queue, clmem_new_p, clmem_p, 0, 0, gaugemomentum_size, NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy old to new gaugemomentum failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::copy_gaugefield_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();

	int clerr = clEnqueueCopyBuffer(queue, clmem_gaugefield, clmem_new_u, 0, 0, sizeof(s_gaugefield), NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy new to old gaugefield failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	(*timer).add();
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	(*timer).reset();
	
	int gaugemomentum_size = sizeof(ae)*GAUGEMOMENTASIZE2;
	int clerr = clEnqueueCopyBuffer(queue, clmem_p, clmem_new_p, 0, 0, gaugemomentum_size, NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy new to old gaugemomentum failed, aborting.";
		exit(HMC_OCLERROR);
	}
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::get_deltah_from_device(hmc_float * out, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	hmc_float tmp;
	copy_float_from_device(clmem_deltah, &tmp, timer);
	(*out) = tmp;
	return HMC_SUCCESS;	
}
