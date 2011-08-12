#include "opencl_hmc.h"
#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

hmc_error Opencl_hmc::fill_collect_options(stringstream* collect_options)
{

	Opencl_fermions::fill_collect_options(collect_options);
	*collect_options <<  " -DBETA=" << get_parameters()->get_beta() << " -DGAUGEMOMENTASIZE=" << GAUGEMOMENTASIZE2;
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fill_buffers()
{
	Opencl_fermions::fill_buffers();

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

	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * SPINORFIELDSIZE;
// 	int eoprec_spinorfield_size = sizeof(spinor) * EOPREC_SPINORFIELDSIZE;
	int gaugemomentum_size = sizeof(ae) * GAUGEMOMENTASIZE2;
// 	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
// 	int global_buf_size = complex_size * num_groups;
// 	int global_buf_size_float = float_size * num_groups;
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
// 	hmc_float tmp;


	/** @todo insert variables needed */
	//init mem-objects

	logger.trace() << "Create buffer for HMC...";
	clmem_force = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);;
	if(clerr != CL_SUCCESS) {
		cout << "creating clmem_force failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clmem_phi_inv = clCreateBuffer(context, CL_MEM_READ_WRITE, spinorfield_size, 0, &clerr);;
	if(clerr != CL_SUCCESS) {
		cout << "creating clmem_phi_inv failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clmem_new_u = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(s_gaugefield), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_new_p = clCreateBuffer(context, CL_MEM_READ_WRITE, gaugemomentum_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_energy_init = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_p2 = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_new_p2 = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_s_fermion = clCreateBuffer(context, CL_MEM_READ_WRITE, float_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}


	return HMC_SUCCESS;
}

void Opencl_hmc::fill_kernels()
{
	//fill kernels of Mother classes
	Opencl_fermions::fill_kernels();

	basic_hmc_code = basic_fermion_code << "types_hmc.h";

	//init kernels for HMC
	set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_hmc_code << "gaugemomentum_zero.cl";
	generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_hmc_code << "random.cl" << "spinorfield_gaussian.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_hmc_code << "random.cl" << "gaugemomentum_gaussian.cl";
	md_update_gaugefield = createKernel("md_update_gaugefield") << basic_hmc_code << "md_update_gaugefield.cl";
	md_update_gaugemomenta = createKernel("md_update_gaugemomenta") << basic_hmc_code << "gaugemomentum.cl" << "md_update_gaugemomenta.cl";
	gauge_force = createKernel("gauge_force") << basic_hmc_code << "gaugemomentum.cl" << "force_gauge.cl";
	fermion_force = createKernel("fermion_force") << basic_hmc_code << "gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_hmc_code << "gaugemomentum_squarenorm.cl";
}

hmc_error Opencl_hmc::init(cl_device_type wanted_device_type, inputparameters* parameters)
{
	hmc_error err = Opencl_fermions::init(wanted_device_type, parameters);
	
	return err;
}

hmc_error Opencl_hmc::finalize_hmc()
{

	logger.debug() << "release HMC-variables.." ;
	if(clReleaseMemObject(clmem_energy_init) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p2) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p2) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_p) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_p) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_new_u) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_phi_inv) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;
	if(clReleaseMemObject(clmem_force) != CL_SUCCESS) return HMC_RELEASEVARIABLEERR;

	logger.debug() << "release HMC-kernels.." ;
	if(clReleaseKernel(generate_gaussian_spinorfield) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	//if(clReleaseKernel(s_gauge)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	//if(clReleaseKernel(s_fermion)!=CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(generate_gaussian_gaugemomenta) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugefield) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(md_update_gaugemomenta) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(gauge_force) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(fermion_force) != CL_SUCCESS) return HMC_RELEASEKERNELERR;
	if(clReleaseKernel(set_zero_gaugemomentum) != CL_SUCCESS) return HMC_RELEASEKERNELERR;

	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Methods needed for the HMC-algorithm

hmc_error Opencl_hmc::generate_gaussian_gaugemomenta_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_new_p
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), &clmem_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( generate_gaussian_gaugemomenta , gs, ls);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::generate_gaussian_spinorfield_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_phi_inv
	clerr = clSetKernelArg(generate_gaussian_spinorfield, 0, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(generate_gaussian_spinorfield, 1, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	enqueueKernel(generate_gaussian_spinorfield  , gs, ls);

	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::md_update_spinorfield_device(const size_t local_work_size, const size_t global_work_size)
{
	//suppose the initial gaussian field is saved in phi_inv. then the "phi" from the algorithm is clmem_inout
	int err =  Opencl_fermions::Qplus_device(clmem_phi_inv, get_clmem_inout() , get_clmem_gaugefield(), local_work_size, global_work_size);

	if(err != HMC_SUCCESS)
		logger.fatal() << "error occured in md_update_spinorfield_device.. ";
	return HMC_SUCCESS;
}


//CP: steps2 is not used at the moment....
hmc_error Opencl_hmc::leapfrog_device(hmc_float tau, int steps1, int steps2, const size_t ls, const size_t gs)
{
	//it is assumed that gaugefield and gaugemomentum have been set to the old ones already
	hmc_float stepsize = (tau) / ((hmc_float) steps1);
	int k;
	hmc_float stepsize_half = 0.5 * stepsize;

	//initial step
	logger.debug() << "\t\tinitial step:";
	//here, phi is inverted using the orig. gaugefield
	force_device(ls, gs);
	md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);
	//intermediate steps
	if(steps1 > 1) logger.debug() << "\t\tperform " << steps1 - 1 << " intermediate steps " ;
	for(k = 1; k < steps1; k++) {
		md_update_gaugefield_device(stepsize, ls, gs);
		force_device(ls, gs);
		md_update_gaugemomentum_device(-1.*stepsize, ls, gs);
	}

	//final step
	logger.debug() << "\t\tfinal step" ;
	md_update_gaugefield_device(stepsize, ls, gs);
	force_device(ls, gs);
	md_update_gaugemomentum_device(-1.*stepsize_half, ls, gs);

	logger.debug() << "\t\tfinished leapfrog";

	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::force_device(const size_t ls, const size_t gs)
{
	logger.debug() << "\t\tstart calculating the force...";
	//CP: make sure that the output field is set to zero
	set_zero_clmem_force_device(ls, gs);
	//add contributions
	logger.debug() << "\t\tcalc gauge_force...";
	gauge_force_device(ls, gs);
// 	cout << "\t\tinvert fermion field..." << endl;
	//CP: to begin with, consider only the cg-solver
	//source is at 0
// 	int k = 0;
// 	int use_cg = TRUE;
	//CP: at the moment, use_eo = 0 so that even-odd is not used!!!!!

	//debugging
//  int err = 0;
//  int use_eo = 0;
//  /** @todo check the use of the sources again, compare to tmlqcd!!! */
//  if(use_cg){
//    if(!use_eo){
//      //the inversion calculates Y = (QplusQminus)^-1 phi = phi_inv
//      cout << "\t\t\tstart solver" << endl;
//      Opencl_fermions::create_point_source_device(k, 0,0,ls, gs, &noop);
//      Opencl_fermions::solver_device(&noop,&noop,&noop,&noop,&noop,&noop,&noop,&noop, ls, gs, get_parameters()->get_cgmax());
//
//      if (err != HMC_SUCCESS) cout << "\t\tsolver did not solve!!" << endl;
//      else cout << "\t\tsolver solved!" << endl;
//    }
//    else{
//      hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
//      hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
//
//      Opencl_fermions::create_point_source_eoprec_device(k,0,0, ls, gs, &noop, &noop, &noop);
//      Opencl_fermions::solver_eoprec_device(&noop,&noop,&noop,&noop,&noop,&noop,&noop,&noop, ls, gs, get_parameters()->get_cgmax());
//    }
//    cout << "\t\t\tcalc X" << endl;
//    //X = Qminus Y = Qminus phi_inv
//    Opencl_fermions::Qminus_device(clmem_phi_inv,get_clmem_inout(), ls, gs, &noop);
//  }
//  else{
//    //here, one has first to invert (with BiCGStab) Qplus phi = X and then invert Qminus X => Qminus^-1 Qplus^-1 phi = (QplusQminus)^-1 phi = Y = phi_inv
//  }
//  cout << "\t\tcalc fermion_force..." << endl;
//  fermion_force_device(ls, gs, &noop);

	return HMC_SUCCESS;
}


hmc_observables Opencl_hmc::metropolis(hmc_float rnd, hmc_float beta, const size_t local_work_size, const size_t global_work_size, usetimer * timer)
{
	//Calc Hamiltonian
	logger.debug() << "Calculate Hamiltonian";
	
	//Gauge-Part
	hmc_float tplaq, splaq, plaq;
	hmc_float tplaq_new, splaq_new, plaq_new;
	hmc_complex poly;
	hmc_complex poly_new;
	//In this call, the observables are calculated already with appropiate Weighting factor of 2.0/(VOL4D*NDIM*(NDIM-1)*NC)
	Opencl::gaugeobservables(get_clmem_gaugefield(), &plaq,  &tplaq, &splaq, &poly);
	Opencl::gaugeobservables(clmem_new_u, &plaq_new,  &tplaq_new, &splaq_new, &poly_new);
	//plaq has to be divided by the norm-factor and multiplied by NC 
	//	(because this is in the defintion of the gauge action and not in the normalization) to get s_gauge
	hmc_float factor = 2.0 / static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) );
	/** NOTE: the minus here is introduced to fit tmlqcd!!! */
	hmc_float deltaH = -(plaq - plaq_new) * beta / factor;
	
	logger.debug() << "\tS_gauge(old field) = " << plaq << "\t" << plaq* beta  /factor;
	logger.debug() << "\tS_gauge(new field) = " << plaq_new << "\t" << plaq_new* beta /factor;
	logger.debug() << "\tdeltaS_gauge = " << deltaH;

	//Gaugemomentum-Part
	hmc_float p2, new_p2;
	set_float_to_gaugemomentum_squarenorm_device(clmem_p, clmem_p2, local_work_size, global_work_size);
	set_float_to_gaugemomentum_squarenorm_device(clmem_new_p, clmem_new_p2, local_work_size, global_work_size);
	Opencl_fermions::copy_float_from_device(clmem_p2, &p2, timer);
	Opencl_fermions::copy_float_from_device(clmem_new_p2, &new_p2, timer);
	//the energy is half the squarenorm
	deltaH += 0.5 * (p2 - new_p2);
	
	logger.debug() << "\tS_gaugemom(old field) = " << 0.5*p2;
	logger.debug() << "\tS_gaugemom(new field) = " << 0.5*new_p2;
	logger.debug() << "\tdeltaS_gaugemom = " << 0.5 * (p2 - new_p2);

	//Fermion-Part:
	hmc_float spinor_energy_init, s_fermion;
	Opencl_fermions::copy_float_from_device(clmem_energy_init, &spinor_energy_init, timer);
	// sum_links phi*_i (M^+M)_ij^-1 phi_j
	// here it is assumed that the rhs has already been computed in clmem_inout... (during the leapfrog..)
	//CP: phi_inv is not needed after this, so it can be used to store M (QplusQminus_inv)
//  Opencl_fermions::Qminus_device(Opencl_fermions::get_clmem_inout(), clmem_phi_inv, local_work_size, global_work_size, timer);
//  set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_s_fermion, local_work_size, global_work_size, timer);
//  copy_float_from_device(clmem_s_fermion, &s_fermion, timer);
//  deltaH += spinor_energy_init - s_fermion;

	//Metropolis-Part
	hmc_float compare_prob;
	if(deltaH < 0) {
		compare_prob = exp(deltaH);
	} else {
		compare_prob = 1.0;
	}
	logger.debug() << "\tdeltaH = " << deltaH << "\tAcc-Prop = " << compare_prob;
	hmc_observables tmp;
	if(rnd <= compare_prob) {
		tmp.accept = 1;
		tmp.plaq = plaq_new;
		tmp.tplaq = tplaq_new;
		tmp.splaq = splaq_new;
		tmp.poly = poly_new;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	} else {
		tmp.accept = 0;
		tmp.plaq = plaq;
		tmp.tplaq = tplaq;
		tmp.splaq = splaq;
		tmp.poly = poly;
		tmp.deltaH = deltaH;
		tmp.prob = compare_prob;
	}

	return tmp;
}

hmc_error Opencl_hmc::calc_spinorfield_init_energy_device(const size_t local_work_size, const size_t global_work_size)
{
	//Suppose the initial spinorfield is saved in phi_inv
	Opencl_fermions::set_float_to_global_squarenorm_device(clmem_phi_inv, clmem_energy_init, local_work_size, global_work_size);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugemomentum_device(hmc_float eps, const size_t ls, const size_t gs)
{
	//__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in){

	hmc_float tmp = eps;

	int clerr;
	clerr = clSetKernelArg(md_update_gaugemomenta, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugemomenta, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugemomenta, 2, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel(md_update_gaugemomenta  , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::md_update_gaugefield_device(hmc_float eps, const size_t ls, const size_t gs)
{
	// __kernel void md_update_gaugefield(hmc_float eps, __global ae * p_in, __global ocl_s_gaugefield * u_inout){
	hmc_float tmp = eps;
	int clerr;
	//this is always applied to clmem_force
	clerr = clSetKernelArg(md_update_gaugefield, 0, sizeof(hmc_float), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugefield, 1, sizeof(cl_mem), &clmem_new_p);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(md_update_gaugefield, 2, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( md_update_gaugefield , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::set_zero_clmem_force_device(const size_t ls, const size_t gs)
{
	int clerr;
	//this is always applied to clmem_force
	clerr = clSetKernelArg(set_zero_gaugemomentum, 0, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( set_zero_gaugemomentum , gs, ls);

	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::gauge_force_device(const size_t ls, const size_t gs)
{
	int clerr;
	clerr = clSetKernelArg(gauge_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(gauge_force, 1, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( gauge_force , gs, ls);
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::fermion_force_device(const size_t ls, const size_t gs)
{
	//fermion_force(field, Y, X, out);
	cl_mem tmp = get_clmem_inout();
	int clerr;
	clerr = clSetKernelArg(fermion_force, 0, sizeof(cl_mem), &clmem_new_u);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 1, sizeof(cl_mem), &clmem_phi_inv);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 2, sizeof(cl_mem), &tmp);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(fermion_force, 3, sizeof(cl_mem), &clmem_force);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel( fermion_force , gs, ls);
	
	return HMC_SUCCESS;

}

hmc_error Opencl_hmc::set_float_to_gaugemomentum_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out, const size_t ls, const size_t gs)
{
	int clerr = CL_SUCCESS;
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
//  /** @todo add reduction */
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), &clmem_out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	
	enqueueKernel(gaugemomentum_squarenorm  , gs, ls);
	
	return HMC_SUCCESS;
}

////////////////////////////////////////////////////
//Methods to copy new and old fields... these can be optimized!!
hmc_error Opencl_hmc::copy_gaugefield_old_new_device(usetimer * timer)
{
	(*timer).reset();
	int gaugefield_size = sizeof(s_gaugefield);
	int clerr = clEnqueueCopyBuffer(queue, clmem_gaugefield, clmem_new_u,  0, 0, gaugefield_size, NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy old to new gaugefield failed, aborting.";
		exit(HMC_OCLERROR);
	}
	(*timer).add();
	
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_old_new_device(usetimer * timer)
{
	(*timer).reset();
	int gaugemomentum_size = sizeof(ae) * GAUGEMOMENTASIZE2;
	int clerr = clEnqueueCopyBuffer(queue, clmem_p, clmem_new_p, 0, 0, gaugemomentum_size, NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy old to new gaugemomentum failed, aborting.";
		exit(HMC_OCLERROR);
	}
	(*timer).add();
	
	return HMC_SUCCESS;
}

hmc_error Opencl_hmc::copy_gaugefield_new_old_device(usetimer * timer)
{
	(*timer).reset();
	int clerr = clEnqueueCopyBuffer(queue, clmem_new_u, clmem_gaugefield, 0, 0, sizeof(s_gaugefield), NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy new to old gaugefield failed, aborting.";
		exit(HMC_OCLERROR);
	}
	(*timer).add();
	
	return HMC_SUCCESS;
}


hmc_error Opencl_hmc::copy_gaugemomenta_new_old_device(usetimer * timer)
{
	(*timer).reset();
	int gaugemomentum_size = sizeof(ae) * GAUGEMOMENTASIZE2;
	int clerr = clEnqueueCopyBuffer(queue, clmem_new_p, clmem_p, 0, 0, gaugemomentum_size, NULL, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy new to old gaugemomentum failed, aborting.";
		exit(HMC_OCLERROR);
	}
	(*timer).add();
	
	return HMC_SUCCESS;
}

#ifdef _PROFILING_
usetimer* Opencl_hmc::get_timer(char * in){
	usetimer *noop = NULL;
	noop = Opencl_fermions::get_timer(in);
	if(noop != NULL) return noop;

	if (strcmp(in, "generate_gaussian_spinorfield") == 0){
    return &this->timer_generate_gaussian_spinorfield;
	}	
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0){
    return &this->timer_generate_gaussian_gaugemomenta;
	}	
	if (strcmp(in, "md_update_gaugefield") == 0){
    return &this->timer_md_update_gaugefield;
	}
	if (strcmp(in, "md_update_gaugemomenta") == 0){
    return &this->timer_md_update_gaugemomenta;
	}
	if (strcmp(in, "gauge_force") == 0){
    return &this->timer_gauge_force;
	}
	if (strcmp(in, "fermion_force") == 0){
    return &this->timer_fermion_force;
	}
	if (strcmp(in, "set_zero_gaugemomentum") == 0){
    return &this->timer_set_zero_gaugemomentum;
	}	
	if (strcmp(in, "gaugemomentum_squarenorm") == 0){
    return &this->timer_gaugemomentum_squarenorm;
	}		
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl_hmc::get_read_write_size(char * in, inputparameters * parameters){
	Opencl::get_read_write_size(in, parameters);
		//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
	  S = EOPREC_SPINORFIELDSIZE;
	else
	  S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "generate_gaussian_spinorfield") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "generate_gaussian_gaugemomenta") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "md_update_gaugefield") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "md_update_gaugemomenta") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "gauge_force") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "fermion_force") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "set_zero_gaugemomentum;") == 0){
    return 10000000000000000000;
	}	
	if (strcmp(in, "gaugemomentum_squarenorm") == 0){
    return 10000000000000000000;
	}	
	
}

void Opencl_hmc::print_profiling(std::string filename){
	Opencl::print_profiling(filename);
	char * kernelName;
	kernelName = "generate_gaussian_spinorfield";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "generate_gaussian_gaugemomenta";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "md_update_gaugefield";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "md_update_gaugemomenta";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gauge_force";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "fermion_force";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_zero_gaugemomentum";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gaugemomentum_squarenorm";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	
}
#endif
