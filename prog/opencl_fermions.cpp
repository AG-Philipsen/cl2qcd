#include "opencl_fermions.h"

#include "logger.hpp"

/**
 * What follows are functions that call opencl_fermions-class-functions.
 * This is needed to be able to pass different fermionmatrices as
 * 	arguments to class-functions.
 */
void M_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf) {
	that->M(in, out, gf);
}
void Qplus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf) {
	that->Qplus(in, out, gf);
}
void Qminus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf) {
	that->Qminus(in, out, gf);
}
void QplusQminus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf) {
	that->QplusQminus(in, out, gf);
}
void Aee_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf) {
	that->Aee(in, out, gf);
}

//Opencl_fermions class functions

hmc_error Opencl_fermions::fill_collect_options(stringstream* collect_options)
{

	Opencl::fill_collect_options(collect_options);
	*collect_options << " -D_FERMIONS_" << " -DSPINORSIZE=" << get_parameters()->get_spinorsize() << " -DHALFSPINORSIZE=" << get_parameters()->get_halfspinorsize() 
		<< " -DSPINORFIELDSIZE=" << get_parameters()->get_spinorfieldsize() << " -DEOPREC_SPINORFIELDSIZE=" << get_parameters()->get_eoprec_spinorfieldsize();
	switch (get_parameters()->get_fermact()) {
		case TWISTEDMASS :
			*collect_options << " -D_TWISTEDMASS_";
			break;
		case CLOVER :
			*collect_options << " -D_CLOVER_";
			break;
	}

	//CP: give kappa and its negative value
	hmc_float kappa_tmp = get_parameters()->get_kappa();
	*collect_options << " -DKAPPA=" << kappa_tmp;
	*collect_options << " -DMKAPPA=" << -kappa_tmp;
	//CP: These are the kappas including BC in spatial and temporal direction
	hmc_float tmp_spatial = (get_parameters()->get_theta_fermion_spatial() * PI) / ( (hmc_float) get_parameters()->get_nspace());
	hmc_float tmp_temporal = (get_parameters()->get_theta_fermion_temporal() * PI) / ( (hmc_float) get_parameters()->get_ntime());
	//BC: on the corners in each direction: exp(i theta) -> on each site exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)
	*collect_options << " -DKAPPA_SPATIAL_RE=" << kappa_tmp*cos(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_RE=" << -kappa_tmp*cos(tmp_spatial);
	*collect_options << " -DKAPPA_SPATIAL_IM=" << kappa_tmp*sin(tmp_spatial);
	*collect_options << " -DMKAPPA_SPATIAL_IM=" << -kappa_tmp*sin(tmp_spatial);

	*collect_options << " -DKAPPA_TEMPORAL_RE=" << kappa_tmp*cos(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_RE=" << -kappa_tmp*cos(tmp_temporal);
	*collect_options << " -DKAPPA_TEMPORAL_IM=" << kappa_tmp*sin(tmp_temporal);
	*collect_options << " -DMKAPPA_TEMPORAL_IM=" << -kappa_tmp*sin(tmp_temporal);

	switch (get_parameters()->get_fermact()) {
		case TWISTEDMASS :
			*collect_options << " -DMU=" << get_parameters()->get_mu();
			get_parameters()->calc_mubar();
			*collect_options << " -DMUBAR=" << get_parameters()->get_mubar();
			get_parameters()->set_mubar_negative();
			*collect_options << " -DMMUBAR=" << get_parameters()->get_mubar();
			get_parameters()->set_mubar_negative();
			break;
		case CLOVER :
			*collect_options << " -DCSW=" << get_parameters()->get_csw();
			break;
	}

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::fill_buffers()
{
	Opencl::fill_buffers();
	
	logger.trace() << "init buffer for solver...";
	int clerr = CL_SUCCESS;

	int spinorfield_size = sizeof(spinor) * SPINORFIELDSIZE;
	int eoprec_spinorfield_size = sizeof(spinor) * EOPREC_SPINORFIELDSIZE;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);

	logger.debug() << "init general spinorfield-buffers";
	clmem_corr = create_rw_buffer(spinorfield_size);
	clmem_inout = create_rw_buffer(spinorfield_size);
	clmem_source = create_rw_buffer(spinorfield_size);
	clmem_tmp = create_rw_buffer(spinorfield_size);
	
	if(get_parameters()->get_use_eo() == true) {
		logger.debug() << "init geberal eoprec-spinorfield-buffers";
		clmem_inout_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_source_even = create_rw_buffer(eoprec_spinorfield_size);
		clmem_source_odd = create_rw_buffer(eoprec_spinorfield_size);
		clmem_tmp_eoprec_1 = create_rw_buffer(eoprec_spinorfield_size);
		clmem_tmp_eoprec_2 = create_rw_buffer(eoprec_spinorfield_size);
		clmem_tmp_eoprec_3 = create_rw_buffer(eoprec_spinorfield_size);
	}
	
	logger.debug() << "init solver spinorfield-buffers";
	///@todo some buffers can be saved here if only cg is used
	if(get_parameters()->get_use_eo() == false) {
		//these are only used in a non-eoprec solver
		clmem_rn = create_rw_buffer(spinorfield_size);
		clmem_rhat = create_rw_buffer(spinorfield_size);
		clmem_v = create_rw_buffer(spinorfield_size);
		clmem_p = create_rw_buffer(spinorfield_size);
		clmem_s =create_rw_buffer(spinorfield_size);
		clmem_t = create_rw_buffer(spinorfield_size);
		clmem_aux = create_rw_buffer(spinorfield_size);
	}
	else{
		//LZ only use the following if we want to apply even odd preconditioning
		logger.debug() << "init solver eoprec-spinorfield-buffers";
		clmem_rn_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_rhat_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_v_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_p_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_s_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_t_eoprec = create_rw_buffer(eoprec_spinorfield_size);
		clmem_aux_eoprec = create_rw_buffer(eoprec_spinorfield_size);

	} //end if: eoprec

	logger.debug() << "create buffers for complex and real numbers";
	clmem_rho = create_rw_buffer(complex_size);
	clmem_rho_next = create_rw_buffer(complex_size);
	clmem_alpha = create_rw_buffer(complex_size);
	clmem_omega = create_rw_buffer(complex_size);
	clmem_beta = create_rw_buffer(complex_size);
	clmem_tmp1 = create_rw_buffer(complex_size);
	clmem_tmp2 = create_rw_buffer(complex_size);
	clmem_one = create_rw_buffer(complex_size);
	clmem_minusone = create_rw_buffer(complex_size);
	clmem_resid = create_rw_buffer(float_size);
	clmem_trueresid = create_rw_buffer(float_size);
	//scratch buffers will be created on demand
	clmem_global_squarenorm_buf_glob = 0;
	clmem_scalar_product_buf_glob = 0;
	
	logger.debug() << "write contents to some buffers";
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	clerr = clEnqueueWriteBuffer(queue, clmem_one, CL_TRUE, 0, complex_size, &one, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... writing clmem_one failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueWriteBuffer(queue, clmem_minusone, CL_TRUE, 0, complex_size, &minusone, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... writing clmem_minusone failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	return HMC_SUCCESS;
}

void Opencl_fermions::fill_kernels()
{
	Opencl::fill_kernels();

	basic_fermion_code = basic_opencl_code << "types_fermions.h" << "operations_su3vec.cl"
	                     << "operations_spinor.cl" << "spinorfield.cl";

	logger.debug() << "Create fermion kernels...";
	if(get_parameters()->get_fermact() == WILSON){
		M_wilson = createKernel("M_wilson") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m.cl";
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		M_tm_plus = createKernel("M_tm_plus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";
		M_tm_minus = createKernel("M_tm_minus") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_m_tm_minus.cl";
	}
	else if(get_parameters()->get_fermact() == CLOVER){
		logger.fatal() << "no kernels for CLOVER-discretization implemented yet, aborting... ";
		exit (HMC_STDERR);
	}
	else{
		logger.fatal() << "there was a problem with which fermion-discretization to use, aborting... ";
		exit (HMC_STDERR);
	}
	
	gamma5 = createKernel("gamma5") << basic_fermion_code << "fermionmatrix.cl" << "fermionmatrix_gamma5.cl";
	
	ps_correlator = createKernel("ps_correlator") << basic_fermion_code << "fermionobservables.cl";

	set_spinorfield_cold = createKernel("set_spinorfield_cold") << basic_fermion_code << "spinorfield_cold.cl";
	saxpy = createKernel("saxpy") << basic_fermion_code << "spinorfield_saxpy.cl";
	saxsbypz = createKernel("saxsbypz") << basic_fermion_code << "spinorfield_saxsbypz.cl";
	scalar_product = createKernel("scalar_product") << basic_fermion_code << "spinorfield_scalar_product.cl";
	scalar_product_reduction = createKernel("scalar_product_reduction") << basic_fermion_code << "spinorfield_scalar_product.cl";
	set_zero_spinorfield = createKernel("set_zero_spinorfield") << basic_fermion_code << "spinorfield_set_zero.cl";
	global_squarenorm = createKernel("global_squarenorm") << basic_fermion_code << "spinorfield_squarenorm.cl";
	global_squarenorm_reduction = createKernel("global_squarenorm_reduction") << basic_fermion_code << "spinorfield_squarenorm.cl";

	ratio = createKernel("ratio") << basic_opencl_code << "complex_ratio.cl";
	product = createKernel("product") << basic_opencl_code << "complex_product.cl";

	convert_to_kappa_format = createKernel("convert_to_kappa_format") << basic_fermion_code << "spinorfield_kappaformat_convert.cl";
	convert_from_kappa_format = createKernel("convert_from_kappa_format") << basic_fermion_code << "spinorfield_kappaformat_convert.cl";
	create_point_source = createKernel("create_point_source") << basic_fermion_code << "spinorfield_point_source.cl";

	//Kernels needed if eoprec is used
	if(get_parameters()->get_use_eo() == true) {
		if(get_parameters()->get_fermact() == TWISTEDMASS){
			M_tm_sitediagonal = createKernel("M_tm_sitediagonal") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
			M_tm_inverse_sitediagonal = createKernel("M_tm_inverse_sitediagonal") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_m.cl";
		}
		dslash_eoprec = createKernel("dslash_eoprec") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo.cl" << "fermionmatrix_eo_dslash.cl";
		gamma5_eoprec = createKernel("gamma5_eoprec") << basic_fermion_code << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "fermionmatrix_eo_gamma5.cl";
		convert_from_eoprec = createKernel("convert_from_eoprec") << basic_fermion_code << "spinorfield_eo_convert.cl";
		set_eoprec_spinorfield_cold = createKernel("set_eoprec_spinorfield_cold") << basic_fermion_code << "spinorfield_eo_cold.cl";
		convert_to_kappa_format_eoprec = createKernel("convert_to_kappa_format_eoprec") << basic_fermion_code << "spinorfield_eo_kappaformat_convert.cl";
		convert_from_kappa_format_eoprec = createKernel("convert_from_kappa_format_eoprec") << basic_fermion_code << "spinorfield_eo_kappaformat_convert.cl";
		saxpy_eoprec = createKernel("saxpy_eoprec") << basic_fermion_code << "spinorfield_eo_saxpy.cl";
		saxsbypz_eoprec = createKernel("saxsbypz_eoprec") << basic_fermion_code << "spinorfield_eo_saxsbypz.cl";
		scalar_product_eoprec = createKernel("scalar_product_eoprec") << basic_fermion_code << "spinorfield_eo_scalar_product.cl";
		set_zero_spinorfield_eoprec = createKernel("set_zero_spinorfield_eoprec") << basic_fermion_code << "spinorfield_eo_zero.cl";
		global_squarenorm_eoprec = createKernel("global_squarenorm_eoprec") << basic_fermion_code << "spinorfield_eo_squarenorm.cl";
		create_point_source_eoprec = createKernel("create_point_source_eoprec") << basic_fermion_code << "spinorfield_eo_point_source.cl";
	}
}

hmc_error Opencl_fermions::init(cl_device_type wanted_device_type, inputparameters* parameters, int nstates)
{
        hmc_error err = Opencl::init(wanted_device_type, parameters, nstates);
	return err;
}

void Opencl_fermions::convert_to_kappa_format_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(convert_to_kappa_format, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_to_kappa_format, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( convert_to_kappa_format, gs2, ls2);
}

void Opencl_fermions::convert_from_kappa_format_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(convert_from_kappa_format, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_kappa_format, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(convert_from_kappa_format, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( convert_from_kappa_format , gs2, ls2);
}

void Opencl_fermions::convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(convert_from_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_eoprec,0,sizeof(cl_mem),&in1); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(convert_from_eoprec,1,sizeof(cl_mem),&in2); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(convert_from_eoprec,2,sizeof(cl_mem),&out); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(convert_from_eoprec , gs2, ls2);
}

	
void Opencl_fermions::convert_to_kappa_format_eoprec_device(cl_mem in){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(convert_to_kappa_format_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_to_kappa_format_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(convert_to_kappa_format_eoprec , gs2, ls2);
}

void Opencl_fermions::convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(convert_from_kappa_format_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_kappa_format_eoprec, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(convert_from_kappa_format_eoprec, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(convert_from_kappa_format_eoprec , gs2, ls2);
}

//compound fermionmatrix-functions without eoprec
void Opencl_fermions::M(cl_mem in, cl_mem out, cl_mem gf){

	if(get_parameters()->get_fermact() == WILSON){
		//in the pure Wilson case there is just one fermionmatrix 
		M_wilson_device(in, out, gf);
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		M_tm_plus_device(in, out, gf);
	}
}

void Opencl_fermions::Qplus(cl_mem in, cl_mem out, cl_mem gf){
	if(get_parameters()->get_fermact() == WILSON){
		//in the pure Wilson case there is just one fermionmatrix 
		M_wilson_device(in, out, gf);
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		M_tm_plus_device(in, out, gf);
	}
	gamma5_device(out);
}

void Opencl_fermions::Qminus(cl_mem in, cl_mem out, cl_mem gf){
	if(get_parameters()->get_fermact() == WILSON){
		//in the pure Wilson case there is just one fermionmatrix 
		M_wilson_device(in, out, gf);
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		M_tm_minus_device(in, out, gf);
	}
	gamma5_device(out);
}

void Opencl_fermions::QplusQminus(cl_mem in, cl_mem out, cl_mem gf)
{
	/** @todo one could save one field here if an additional copying would be included in the end... 
	 * or the field should be created in here, local */
	Qminus(in, clmem_tmp, gf);
	Qplus(clmem_tmp, out, gf);
}

//explicit fermionmatrix-kernel calling functions
void Opencl_fermions::M_wilson_device(cl_mem in, cl_mem out, cl_mem gf){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(M_wilson, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
  int clerr = clSetKernelArg(M_wilson,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_wilson,1,sizeof(cl_mem),&gf);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_wilson,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel( M_wilson, gs2, ls2);
}

void Opencl_fermions::M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(M_tm_plus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
  int clerr = clSetKernelArg(M_tm_plus,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_tm_plus,1,sizeof(cl_mem),&gf);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_tm_plus,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel( M_tm_plus, gs2, ls2);
}

void Opencl_fermions::M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(M_tm_minus, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
  int clerr = clSetKernelArg(M_tm_minus,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_tm_minus,1,sizeof(cl_mem),&gf);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_tm_minus,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel( M_tm_minus, gs2, ls2);
}

void Opencl_fermions::gamma5_device(cl_mem inout){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(gamma5, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(gamma5 , gs2, ls2);
}

//compound fermionmatrix-functions with eoprec
void Opencl_fermions::Aee(cl_mem in, cl_mem out, cl_mem gf)
{
	int even = EVEN;
	int odd = ODD;

	if(get_parameters()->get_fermact() == WILSON){
		logger.fatal() << "not yet implemented in WILSON case, aborting...";
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		dslash_eoprec_device(in, clmem_tmp_eoprec_1, gf, odd);
		M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2);
		dslash_eoprec_device(clmem_tmp_eoprec_2, out, gf, even);
		M_tm_sitediagonal_device(in, clmem_tmp_eoprec_1);
		/** @todo  the copying can be extincted by using
		 * 	saxpy_eoprec_device(out, clmem_tmp_eoprec_1, clmem_one, out);
		 */
		copy_buffer_on_device(out, clmem_tmp_eoprec_3, sizeof(spinor) * EOPREC_SPINORFIELDSIZE);
		saxpy_eoprec_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1, clmem_one, out);
	}
}

//explicit eoprec fermionmatrix functions
void Opencl_fermions::gamma5_eoprec_device(cl_mem inout){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(gamma5_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(gamma5_eoprec,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel( gamma5_eoprec, gs2, ls2);
}

void Opencl_fermions::dslash_eoprec_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd)
{
	int eo = evenodd;
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(dslash_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(dslash_eoprec,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,2,sizeof(cl_mem),&gf);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash_eoprec,3,sizeof(int),&eo);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	enqueueKernel(dslash_eoprec , gs2, ls2);
}

void Opencl_fermions::M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(M_tm_inverse_sitediagonal, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(M_tm_inverse_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( M_tm_inverse_sitediagonal, gs2, ls2);
}

void Opencl_fermions::M_tm_sitediagonal_device(cl_mem in, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(M_tm_sitediagonal, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(M_tm_sitediagonal, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(M_tm_sitediagonal, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(M_tm_sitediagonal , gs2, ls2);
}

//BLAS-functions
void Opencl_fermions::saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(saxpy, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(saxpy , gs2, ls2);
}

void Opencl_fermions::set_spinorfield_cold_device(cl_mem inout){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_spinorfield_cold, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_spinorfield_cold,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(set_spinorfield_cold , gs2, ls2);
}

void Opencl_fermions::set_eoprec_spinorfield_cold_device(cl_mem inout){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_eoprec_spinorfield_cold, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_eoprec_spinorfield_cold,0,sizeof(cl_mem),&inout); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(set_eoprec_spinorfield_cold , gs2, ls2);
}

void Opencl_fermions::saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(saxpy_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxpy_eoprec, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( saxpy_eoprec, gs2, ls2);
}

void Opencl_fermions::saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(saxsbypz, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( saxsbypz, gs2, ls2);
}

void Opencl_fermions::saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(saxsbypz_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(saxsbypz_eoprec, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( saxsbypz_eoprec, gs2, ls2);
}

void Opencl_fermions::set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(scalar_product, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_complex = sizeof(hmc_complex)*num_groups;
	if( clmem_scalar_product_buf_glob == 0 ) clmem_scalar_product_buf_glob = create_rw_buffer(global_buf_size_complex);

	//set arguments
	int clerr = clSetKernelArg(scalar_product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}

	clerr = clSetKernelArg(scalar_product, 3, sizeof(hmc_complex) * ls2, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(scalar_product , gs2, ls2);

	/** @todo Here the wait is needed. Replace this call by a clWaitForEvents! */
	clFinish(queue);
	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( scalar_product_reduction, gs2, ls2);
}

void Opencl_fermions::set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(scalar_product_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_complex = sizeof(hmc_complex)*num_groups;
	if( clmem_scalar_product_buf_glob == 0 ) clmem_scalar_product_buf_glob = create_rw_buffer(global_buf_size_complex);	
	
	//set arguments
	int clerr = clSetKernelArg(scalar_product_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_eoprec, 3, sizeof(hmc_complex) * ls2, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( scalar_product_eoprec, gs2, ls2);
	clFinish(queue);
	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( scalar_product_reduction, gs2, ls2);
}


void Opencl_fermions::set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(ratio, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ratio, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(ratio, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(ratio, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP:this needs only one kernel!!
	ls2 = 1;
	gs2 = 1;
	enqueueKernel( ratio, gs2, ls2);
}

void Opencl_fermions::set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(product, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(product, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP:this needs only one kernel!!
	ls2 = 1;
	gs2 = 1;
	enqueueKernel(product , gs2, ls2);
}

void Opencl_fermions::set_float_to_global_squarenorm_device(cl_mem a, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(global_squarenorm, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_float = sizeof(hmc_float)*num_groups;
	if( clmem_global_squarenorm_buf_glob == 0 ) clmem_global_squarenorm_buf_glob = create_rw_buffer(global_buf_size_float);	
	
	//set arguments
	int clerr = clSetKernelArg(global_squarenorm, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm, 2, sizeof(hmc_float) * ls2, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(global_squarenorm , gs2, ls2);
	clFinish(queue);
	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( global_squarenorm_reduction, gs2, ls2);
}

void Opencl_fermions::set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(global_squarenorm_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int global_buf_size_float = sizeof(hmc_float)*num_groups;
	if( clmem_global_squarenorm_buf_glob == 0 ) clmem_global_squarenorm_buf_glob = create_rw_buffer(global_buf_size_float);	
	
	int clerr = clSetKernelArg(global_squarenorm_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_eoprec, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_eoprec, 2, sizeof(hmc_float) * ls2, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( global_squarenorm_eoprec, gs2, ls2);
	clFinish(queue);
	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( global_squarenorm_reduction, gs2, ls2);
}

void Opencl_fermions::set_zero_spinorfield_device(cl_mem x)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_zero_spinorfield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(set_zero_spinorfield , gs2, ls2);
}

void Opencl_fermions::set_zero_spinorfield_eoprec_device(cl_mem x)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_zero_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( set_zero_spinorfield_eoprec, gs2, ls2);
}

hmc_error Opencl_fermions::bicgstab_device( matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, int cgmax)
{

int debug = 0;
if(debug) cout << "debug-output at bicgstab_device is activated" << endl;
	
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter = 0; iter < cgmax; iter++) {
		if(iter % iter_refresh == 0) {
			set_zero_spinorfield_device(clmem_v);
			set_zero_spinorfield_device(clmem_p);
			f(this, inout, clmem_rn, gf);

			saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
			copy_buffer_on_device(clmem_rn, clmem_rhat, sizeof(spinor) * SPINORFIELDSIZE);

			copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
			copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
			copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));

			//CP: calc initial residuum for output, this is not needed for the algorithm!!
// 			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size);
//			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
// 			cout << "initial residuum at iter " << iter << " is: " << scientific << resid << endl;
			//printf("initial residuum at iter %i is %.40e\n", iter, resid);
		}

////////////////////////////////////
//collect all variables if debug is enabled
if(debug){
	hmc_complex omega;
	hmc_complex rho;
	hmc_complex rho_next;
	hmc_complex tmp1;
	hmc_complex tmp2;
	hmc_complex beta;
	hmc_complex alpha;
		
	get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
	get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
	get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));
	
 	cout << "debug output at start: " << endl;
	cout << " rho: " << rho.re << "  " <<  rho.im << endl;
	cout << " alpha: " <<alpha.re << "  " <<  alpha.im << endl;
	cout << " omega: " << omega.re << "  " <<  omega.im << endl;
}
////////////////////////////////////

		set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next);

if(debug){
	hmc_complex rho_next;
	get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
		cout << "(rhat, rn): " << rho_next.re << "  " <<  rho_next.im << endl;
}
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
		copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
		saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p);

if(debug){
	hmc_complex rho_next;
	set_complex_to_scalar_product_device(clmem_p, clmem_p, clmem_tmp2);
	get_buffer_from_device(clmem_tmp2, &rho_next, sizeof(hmc_complex));
		cout << "(p,p): " << rho_next.re << "  " <<  rho_next.im << endl;
				set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
}

		f(this, clmem_p, clmem_v, gf);

////////////////////////////////////
//collect all variables if debug is enabled
if(debug){
	hmc_complex omega;
	hmc_complex rho;
	hmc_complex rho_next;
	hmc_complex tmp1;
	hmc_complex tmp2;
	hmc_complex beta;
	hmc_complex alpha;
		
	get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
	get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
	get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));
	
 	cout << "debug output after first half: " << endl;
	cout << " rho: " << rho.re << "  " <<  rho.im << endl;
	cout << " rho_next: " << rho_next.re << "  " <<  rho_next.im << endl;
	cout << " beta: " << beta.re << "  " <<  beta.im << endl;
	cout << " alpha: " <<alpha.re << "  " <<  alpha.im << endl;
	cout << " omega: " << omega.re << "  " <<  omega.im << endl;
	cout << " tmp1: " << tmp1.re << "  " <<  tmp1.im << endl;
	cout << " tmp2: " << tmp2.re << "  " <<  tmp2.im << endl;
}
////////////////////////////////////


		set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1);

if(debug){
	hmc_complex rho_next;
	get_buffer_from_device(clmem_tmp1, &rho_next, sizeof(hmc_complex));
		cout << "(rhat, v): " << rho_next.re << "  " <<  rho_next.im << endl;
}

		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
		saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s);

// 		//see if s is too small
// 		hmc_complex s_norm;
// 		//borrow clmem_alpha for this
// 		set_complex_to_scalar_product_device(clmem_s, clmem_s, clmem_alpha, localsize, globalsize);
// 		get_buffer_from_device(clmem_alpha, &s_norm, sizeof(hmc_complex));
// 		if(debug) cout << "|s|^2: " << s_norm.re << "  " <<  s_norm.im << endl;
// 		//reset value of alpha
// 		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
// 		//check if |s|^2 is too small
// 		if(s_norm.re < epssquare){
// 			set_complex_to_product_device(clmem_minusone, clmem_alpha, clmem_alpha);
// 			saxpy_device(clmem_p, inout, clmem_alpha, inout, localsize, globalsize);
// 			
// 			f(this, inout, clmem_aux, gf, localsize, globalsize);
// 			saxpy_device(clmem_aux, source, clmem_one, clmem_aux, localsize, globalsize);
// 			set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid, localsize, globalsize);
// 			get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
// 			cout << "\ttrueresiduum:\t" << trueresid << " has to be smaller than " << epssquare << endl;
// 
// 			cout << "|s|^2 is too small to continue..." << endl;
// 
// 			return HMC_SUCCESS;
// 		}
		

		f(this, clmem_s, clmem_t, gf);

		set_complex_to_scalar_product_device(clmem_t, clmem_s, clmem_tmp1);
		//!!CP: this can also be global_squarenorm
		set_complex_to_scalar_product_device(clmem_t, clmem_t, clmem_tmp2);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);

		saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn);

		saxsbypz_device(clmem_p, clmem_s, inout, clmem_alpha, clmem_omega, inout);

if(debug){
	hmc_complex rho_next;
	set_complex_to_scalar_product_device(inout, inout, clmem_omega);
	get_buffer_from_device(clmem_omega, &rho_next, sizeof(hmc_complex));
		cout << "(inout,inout): " << rho_next.re << "  " <<  rho_next.im << endl;
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);
}				
		
////////////////////////////////////
//collect all variables if debug is enabled
if(debug){
	hmc_complex omega;
	hmc_complex rho;
	hmc_complex rho_next;
	hmc_complex tmp1;
	hmc_complex tmp2;
	hmc_complex beta;
	hmc_complex alpha;
		
	get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
	get_buffer_from_device(clmem_rho_next, &rho_next, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
	get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
	get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
	get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));
	
 	cout << "debug output after second half: " << endl;
	cout << " rho: " << rho.re << "  " <<  rho.im << endl;
	cout << " rho_next: " << rho_next.re << "  " <<  rho_next.im << endl;
	cout << " beta: " << beta.re << "  " <<  beta.im << endl;
	cout << " alpha: " <<alpha.re << "  " <<  alpha.im << endl;
	cout << " omega: " << omega.re << "  " <<  omega.im << endl;

	cout << " tmp1: " << tmp1.re << "  " <<  tmp1.im << endl;
	cout << " tmp2: " << tmp2.re << "  " <<  tmp2.im << endl;
}
////////////////////////////////////

		set_float_to_global_squarenorm_device(clmem_rn, clmem_resid);
		get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		cout << "resid at iter " << iter << " is: " << resid << endl;
		
		if(resid < epssquare) {
			f(this, inout, clmem_aux, gf);
			saxpy_device(clmem_aux, source, clmem_one, clmem_aux);
			set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid);
			get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
			cout << "\tsolver converged! residuum:\t" << resid << " is smaller than " << epssquare << endl;
			cout << "\ttrueresiduum:\t" << trueresid << " has to be smaller than " << epssquare << endl;
			if(trueresid < epssquare)
				return HMC_SUCCESS;
			else {
				cout << "trueresiduum not small enough" <<endl;
// 				hmc_complex s_norm;
// 				//borrow clmem_alpha for this
// 				set_complex_to_scalar_product_device(clmem_s, clmem_s, clmem_alpha, localsize, globalsize);
// 				get_buffer_from_device(clmem_alpha, &s_norm, sizeof(hmc_complex));
// 				cout << "|s|^2: " << s_norm.re << "  " <<  s_norm.im << endl;
// 				//reset value of alpha
// 				set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);
// 				//check if |s|^2 is too small
// 				if(s_norm.re < epssquare){
// 					cout << "|s|^2 is too small to continue..." << endl;
// // 					return HMC_SUCCESS;
// 				}
			}
		} 
		else {
			printf("residuum at iter%i is:\t%.10e\n", iter, resid);//cout << "residuum:\t" << resid << endl;
		}
	}
	return HMC_STDERR;
}

hmc_error Opencl_fermions::bicgstab_eoprec_device(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, int cgmax)
{
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	for(int iter = 0; iter < cgmax; iter++) {

		if(iter % iter_refresh == 0) {
			set_zero_spinorfield_eoprec_device(clmem_v_eoprec);
			set_zero_spinorfield_eoprec_device(clmem_p_eoprec);

			f(this, inout, clmem_rn_eoprec, gf);

			saxpy_eoprec_device(clmem_rn_eoprec, source, clmem_one, clmem_rn_eoprec);
			copy_buffer_on_device(clmem_rn_eoprec, clmem_rhat_eoprec, sizeof(spinor) * EOPREC_SPINORFIELDSIZE);
			
			copy_buffer_on_device(clmem_one, clmem_alpha, sizeof(hmc_complex));
			copy_buffer_on_device(clmem_one, clmem_omega, sizeof(hmc_complex));
			copy_buffer_on_device(clmem_one, clmem_rho, sizeof(hmc_complex));

			//!!CP: calc initial residuum, this is not needed for the algorithm!!
			//set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
			//get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
			//cout << "initial residuum is: " << resid << endl;
		}

		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1);
		copy_buffer_on_device(clmem_rho_next, clmem_rho, sizeof(hmc_complex));
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2);
		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec);

		f(this, clmem_p_eoprec, clmem_v_eoprec, gf);

		set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_v_eoprec, clmem_tmp1);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha);

		saxpy_eoprec_device(clmem_v_eoprec, clmem_rn_eoprec, clmem_alpha, clmem_s_eoprec);

		f(this, clmem_s_eoprec, clmem_t_eoprec, gf);

		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_tmp1);
		//!!CP: can this also be global_squarenorm??
		set_complex_to_scalar_product_eoprec_device(clmem_t_eoprec, clmem_t_eoprec, clmem_tmp2);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega);

		saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec);

		saxsbypz_eoprec_device(clmem_p_eoprec, clmem_s_eoprec, inout, clmem_alpha, clmem_omega, inout);

		set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid);
		get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));

		if(resid < epssquare) {
			f(this, inout, clmem_aux_eoprec, gf);
			saxpy_eoprec_device(clmem_aux_eoprec, clmem_source_even, clmem_one, clmem_aux_eoprec);
			set_float_to_global_squarenorm_eoprec_device(clmem_aux_eoprec, clmem_trueresid);
			get_buffer_from_device(clmem_trueresid, &trueresid, sizeof(hmc_float));
			//cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid < epssquare)
				return HMC_SUCCESS;
		} else {
			//      cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::cg_device(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, int cgmax)
{
	//CP: these have to be on the host
	hmc_float resid;
	int iter;
	for(iter = 0; iter < cgmax; iter ++) {
		if(iter % iter_refresh == 0) {
			f(this, inout, clmem_rn, gf);
			saxpy_device(clmem_rn, source, clmem_one, clmem_rn);
			copy_buffer_on_device(clmem_rn, clmem_p, sizeof(spinor) * SPINORFIELDSIZE);

		}
		//alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_device(clmem_rn, clmem_rn, clmem_omega);
// 		hmc_complex omega;
// 		get_buffer_from_device(clmem_omega, &omega, sizeof(hmc_complex));
// 		cout << "omega: " << omega.re << " " << omega.im << endl;
		//A pn --> v
		f(this,clmem_p, clmem_v, gf);
		set_complex_to_scalar_product_device(clmem_p, clmem_v, clmem_rho);
// 		hmc_complex rho;
// 		get_buffer_from_device(clmem_rho, &rho, sizeof(hmc_complex));
// 		cout << "rho: " << rho.re << " " << rho.im << endl;
		
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha);	
// 		hmc_complex alpha;
// 		get_buffer_from_device(clmem_alpha, &alpha, sizeof(hmc_complex));
// 		cout << "alpha: " << alpha.re << " " << alpha.im << endl;
		
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1);
// 		hmc_complex tmp1;
// 		get_buffer_from_device(clmem_tmp1, &tmp1, sizeof(hmc_complex));
// 		cout << "tmp1: " << tmp1.re << " " << tmp1.im << endl;
		
		//xn+1
		saxpy_device(inout, clmem_p, clmem_tmp1, inout);
		//rn+1 -> rhat
		saxpy_device(clmem_rn, clmem_v, clmem_alpha, clmem_rhat);

		set_float_to_global_squarenorm_device(clmem_rhat, clmem_resid);
		get_buffer_from_device(clmem_resid, &resid, sizeof(hmc_float));
		cout << "resid: " << resid << endl;

		if(resid < epssquare) {
			//???
			//copy_buffer_on_device(clmem_rhat, clmem_inout, sizeof(spinor) * SPINORFIELDSIZE);

			return HMC_SUCCESS;
		} else {
			//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_rho_next);
			set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta);
// 					hmc_complex beta;
// 		get_buffer_from_device(clmem_beta, &beta, sizeof(hmc_complex));
// 		cout << "beta: " << beta.re << " " << beta.im << endl;

			//pn+1 = rn+1 + beta*pn
			set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2);
// 		hmc_complex tmp2;
// 		get_buffer_from_device(clmem_tmp2, &tmp2, sizeof(hmc_complex));
// 		cout << "tmp2: " << tmp2.re << " " << tmp2.im << endl;
			
			saxpy_device(clmem_p, clmem_rhat, clmem_tmp2, clmem_p);
			
			//rn = rn+1 ^= rn = rhat
			copy_buffer_on_device(clmem_rhat, clmem_rn, sizeof(spinor) * SPINORFIELDSIZE);

		}
	}
	return HMC_SUCCESS;
}

hmc_error Opencl_fermions::cg_eoprec_device(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, int cgmax){
	/// to be implemented if the above one has been checked..
	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::solver_device(matrix_function_call f, cl_mem inout, cl_mem source, cl_mem gf, usetimer * solvertimer, int cgmax)
{
	(*solvertimer).reset();
	convert_to_kappa_format_device(clmem_inout);
	convert_to_kappa_format_device(clmem_source);
	if(get_parameters()->get_use_cg() == true)
	 	cg_device(f, inout, source, gf, cgmax);
	else 
		bicgstab_device(f, inout, source, gf, cgmax);
	convert_from_kappa_format_device(clmem_inout, clmem_inout);
	convert_from_kappa_format_device(clmem_source, clmem_source);
	clFinish(queue);
	(*solvertimer).add();

	return HMC_SUCCESS;
}


hmc_error Opencl_fermions::solver_eoprec_device(matrix_function_call f, cl_mem inout, cl_mem inout_eo, cl_mem source_even, cl_mem source_odd, cl_mem gf, usetimer * solvertimer, int cgmax)
{
	(*solvertimer).reset();

	//CP: even solution
	convert_to_kappa_format_eoprec_device(clmem_inout_eoprec);
	convert_to_kappa_format_eoprec_device(clmem_source_even);
	if(get_parameters()->get_use_cg() == true)
	 	cg_eoprec_device(f, inout_eo, source_even, gf, cgmax);
	else 
		bicgstab_eoprec_device(f, inout_eo, source_even, gf, cgmax);

	//P: odd solution
	/** @todo CP: perhaps one can save some variables used here */
	if(get_parameters()->get_fermact() == WILSON){
		logger.fatal() << "not yet implemented in WILSON case, aborting...";
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		///@todo here, the use of clmem_tmp_eoprec_3 can be avoided
		dslash_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_3, gf, ODD);
		M_tm_inverse_sitediagonal_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_1);
		M_tm_inverse_sitediagonal_device(clmem_source_odd, clmem_tmp_eoprec_2);
		saxpy_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_2, clmem_one, clmem_tmp_eoprec_3);

		convert_from_kappa_format_eoprec_device(clmem_tmp_eoprec_3, clmem_tmp_eoprec_3);
		convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec);
		//CP: whole solution
		//CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eoprec_3
		convert_from_eoprec_device(clmem_inout_eoprec, clmem_tmp_eoprec_3, clmem_inout);
		///convert source back from kappa???
	}
	clFinish(queue);
	(*solvertimer).add();
	return HMC_SUCCESS;
}

void Opencl_fermions::create_point_source_device(cl_mem inout, int i, int spacepos, int timepos)
{
	set_zero_spinorfield_device(inout);

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_zero_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg 3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel( create_point_source, gs2, ls2);
}


void Opencl_fermions::create_point_source_eoprec_device(cl_mem inout_even, cl_mem inout_odd, cl_mem gf, int i, int spacepos, int timepos)
{
	set_zero_spinorfield_eoprec_device(inout_even);
	set_zero_spinorfield_eoprec_device(inout_odd);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_2);
	set_zero_spinorfield_eoprec_device(clmem_tmp_eoprec_1);

	int glob_pos = get_global_pos(spacepos, timepos);
	int n = get_n_eoprec(spacepos, timepos);
	int clerr = CL_SUCCESS;
	int evenodd = glob_pos % 2;

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_zero_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//CP: this is different than the host code, where this is done implicitly when converting the normal source to even/odd
	if(evenodd == 1) {
		clerr = clSetKernelArg(create_point_source_eoprec, 0, sizeof(cl_mem), &inout_odd);
		if(clerr != CL_SUCCESS) {
			cout << "clSetKernelArg 0 failed with " << clerr << " , aborting..." << endl;
			exit(HMC_OCLERROR);
		}
  }
	else {
		clerr = clSetKernelArg(create_point_source_eoprec,0,sizeof(cl_mem),&clmem_tmp_eoprec_2);
		if(clerr!=CL_SUCCESS) {
		  cout<<"clSetKernelArg 0 failed with " << clerr << " , aborting..."<<endl;
    exit(HMC_OCLERROR);
		}		
	}
	clerr = clSetKernelArg(create_point_source_eoprec,1,sizeof(int),&i);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(create_point_source_eoprec,2,sizeof(int),&n);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	enqueueKernel(create_point_source_eoprec , gs2, ls2);

	if(get_parameters()->get_fermact() == WILSON){
		logger.fatal() << "not yet implemented in WILSON case, aborting...";
	}
	else if(get_parameters()->get_fermact() == TWISTEDMASS){
		///@todo here, the use of clmem_tmp_eoprec_3 can be avoided by using clmem_source_even instead
		M_tm_inverse_sitediagonal_device(inout_odd, clmem_tmp_eoprec_1);
		dslash_eoprec_device(clmem_tmp_eoprec_1, clmem_tmp_eoprec_3, gf, EVEN);
		saxpy_eoprec_device(clmem_tmp_eoprec_2, clmem_tmp_eoprec_3, clmem_one, inout_even);
	}
}

//functions to calculate the correlator
void Opencl_fermions::ps_correlator_device(cl_mem in){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes2(set_zero_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ps_correlator,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  ///@todo improve this
	ls2 = 1;
	gs2 = 1;
	enqueueKernel(ps_correlator , gs2, ls2);
}

cl_mem Opencl_fermions::get_clmem_inout()
{
	return clmem_inout;
}

cl_mem Opencl_fermions::get_clmem_source()
{
	return clmem_source;
}

cl_mem Opencl_fermions::get_clmem_tmp()
{
	return clmem_tmp;
}

cl_mem Opencl_fermions::get_clmem_corr()
{
	return clmem_corr;
}

cl_mem Opencl_fermions::get_clmem_inout_eoprec()
{
	return clmem_inout_eoprec;
}

cl_mem Opencl_fermions::get_clmem_tmp_eoprec_1()
{
	return clmem_tmp_eoprec_1;
}

cl_mem Opencl_fermions::get_clmem_source_even()
{
	return clmem_source_even;
}

cl_mem Opencl_fermions::get_clmem_source_odd()
{
	return clmem_source_odd;
}

cl_mem Opencl_fermions::get_clmem_minusone()
{
	return clmem_minusone;
}

hmc_error Opencl_fermions::finalize_fermions()
{

	if(clReleaseKernel(ps_correlator) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(set_spinorfield_cold) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(M_wilson) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(saxpy) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(saxsbypz) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(scalar_product) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(scalar_product_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);

	if(clReleaseKernel(set_zero_spinorfield) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(global_squarenorm) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(global_squarenorm_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(ratio) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(product) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_to_kappa_format) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(convert_from_kappa_format) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(create_point_source) != CL_SUCCESS) exit(HMC_OCLERROR);

	if(get_parameters()->get_use_eo() == false) {
		if(clReleaseKernel(convert_to_kappa_format_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(convert_from_kappa_format_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(dslash_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(saxpy_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(scalar_product_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(set_zero_spinorfield_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(global_squarenorm_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(create_point_source_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(M_tm_sitediagonal) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(M_tm_inverse_sitediagonal) != CL_SUCCESS) exit(HMC_OCLERROR);
	}

	if(clReleaseMemObject(clmem_inout) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_source) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rn) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rhat) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_v) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_p) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_s) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_t) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_aux) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp) != CL_SUCCESS) exit(HMC_OCLERROR);

	if(get_parameters()->get_use_eo() == false) {
		if(clReleaseMemObject(clmem_inout_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_source_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_source_odd) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_rn_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_rhat_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_v_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_p_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_s_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_t_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_aux_eoprec) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_tmp_eoprec_1) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_tmp_eoprec_2) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseMemObject(clmem_tmp_eoprec_3) != CL_SUCCESS) exit(HMC_OCLERROR);
	}

	if(clReleaseMemObject(clmem_rho) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rho_next) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_alpha) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_omega) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_beta) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp1) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tmp2) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_one) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_minusone) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_scalar_product_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_global_squarenorm_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_resid) != CL_SUCCESS) exit(HMC_OCLERROR);
	/** @todo CP: this call seems to cause a seg fault!! */
	if(clReleaseMemObject(clmem_trueresid) != CL_SUCCESS) exit(HMC_OCLERROR);

	return HMC_SUCCESS;
}

#ifdef _PROFILING_
usetimer* Opencl_fermions::get_timer(char * in){
	usetimer *noop = NULL;
	noop = Opencl::get_timer(in);
	if(noop != NULL) return noop;
	
	if (strcmp(in, "M_wilson") == 0){
    return &(this->timer_M_wilson);
	}	
	if (strcmp(in, "gamma5") == 0){
    return &this->timer_gamma5;
	}
	if (strcmp(in, "M_tm_plus") == 0){
    return &this->timer_M_tm_plus;
	}
	if (strcmp(in, "M_tm_minus") == 0){
    return &this->timer_M_tm_minus;
	}
	if (strcmp(in, "gamma5_eoprec") == 0){
    return &this->timer_gamma5_eoprec;
	}
	if (strcmp(in, "dslash_eoprec") == 0){
    return &this->timer_dslash_eoprec;
	}
	if (strcmp(in, "set_spinorfield_cold") == 0){
    return &this->timer_set_spinorfield_cold;
	}
	if (strcmp(in, "set_eoprec_spinorfield_cold") == 0){
    return &this->timer_set_eoprec_spinorfield_cold;
	}
	if (strcmp(in, "convert_from_eoprec") == 0){
    return &this->timer_convert_from_eoprec;
	}
	if (strcmp(in, "saxpy") == 0){
    return &(this->timer_saxpy);
	}
	if (strcmp(in, "saxsbypz") == 0){
    return &this->timer_saxsbypz;
	}
	if (strcmp(in, "set_zero_spinorfield") == 0){
    return &this->timer_set_zero_spinorfield;
	}
	if (strcmp(in, "convert_to_kappa_format") == 0){
    return &this->timer_convert_to_kappa_format;
	}
	if (strcmp(in, "convert_from_kappa_format") == 0){
    return &this->timer_convert_from_kappa_format;
	}
	if (strcmp(in, "convert_to_kappa_format_eoprec") == 0){
    return &this->timer_convert_to_kappa_format_eoprec;
	}
	if (strcmp(in, "convert_from_kappa_format_eoprec") == 0){
    return &this->timer_convert_from_kappa_format_eoprec;
	}
	if (strcmp(in, "create_point_source") == 0){
    return &this->timer_create_point_source;
	}
	if (strcmp(in, "saxpy_eoprec") == 0){
    return &this->timer_saxpy_eoprec;
	}
	if (strcmp(in, "saxsbypz_eoprec") == 0){
    return &this->timer_saxsbypz_eoprec;
	}
	if (strcmp(in, "set_zero_spinorfield_eoprec") == 0){
    return &this->timer_set_zero_spinorfield_eoprec;
	}
	if (strcmp(in, "create_point_source_eoprec") == 0){
    return &this->timer_create_point_source_eoprec;
	}
	if (strcmp(in, "scalar_product") == 0){
    return &this->timer_scalar_product;
	}
	if (strcmp(in, "scalar_product_reduction") == 0){
    return &this->timer_scalar_product_reduction;
	}
	if (strcmp(in, "global_squarenorm") == 0){
    return &this->timer_global_squarenorm;
	}
	if (strcmp(in, "global_squarenorm_reduction") == 0){
    return &this->timer_global_squarenorm_reduction;
	}
	if (strcmp(in, "scalar_product_eoprec") == 0){
    return &this->timer_scalar_product_eoprec;
	}
	if (strcmp(in, "global_squarenorm_eoprec") == 0){
    return &this->timer_global_squarenorm_eoprec;
	}
	if (strcmp(in, "ratio") == 0){
    return &this->timer_ratio;
	}
	if (strcmp(in, "product") == 0){
    return &this->timer_product;
	}
	if (strcmp(in, "ps_correlator") == 0){
    return &this->timer_ps_correlator;
	}
	
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl_fermions::get_read_write_size(char * in, inputparameters * parameters){
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
	if (strcmp(in, "M_wilson") == 0){
    return (240 + 16*R)*D*S;
	}	
	if (strcmp(in, "gamma5") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "M_tm_plus") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "M_tm_minus") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "gamma5_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "M_tm_sitediagonal") == 0){
    return 48*D*S;
	}
	if (strcmp(in, "M_tm_inverse_sitediagonal") == 0){
    return 48*D*S;
	}
	if (strcmp(in, "dslash_eoprec") == 0){
    return (216 + 16*R)*D*S;
	}
	if (strcmp(in, "set_spinorfield_cold") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "set_eoprec_spinorfield_cold") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "convert_from_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "saxpy") == 0){
    return 74*D*S;
	}
	if (strcmp(in, "saxsbypz") == 0){
    return 100*D*S;
	}
	if (strcmp(in, "set_zero_spinorfield") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "convert_to_kappa_format") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "convert_from_kappa_format") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "convert_to_kappa_format_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "convert_from_kappa_format_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "create_point_source") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "saxpy_eoprec") == 0){
    return 74*D*S;
	}
	if (strcmp(in, "saxsbypz_eoprec") == 0){
    return 100*D*S;
	}
	if (strcmp(in, "set_zero_spinorfield_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "create_point_source_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "scalar_product") == 0){
    return 50*D*S;
	}
	if (strcmp(in, "scalar_product_reduction") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "global_squarenorm") == 0){
    return 25*D*S;
	}
	if (strcmp(in, "global_squarenorm_reduction") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "scalar_product_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "global_squarenorm_eoprec") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "ratio") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "product") == 0){
    return 1000000000000000000000000;
	}
	if (strcmp(in, "ps_correlator") == 0){
    return 1000000000000000000000000;
	}
	return 0;
}

void Opencl_fermions::print_profiling(std::string filename){
	Opencl::print_profiling(filename);
	char * kernelName;
	kernelName = "M_wilson";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gamma5";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "M_tm_plus";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "M_tm_minus";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "gamma5_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "M_tm_sitediagonal";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "M_tm_inverse_sitediagonal";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "dslash_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_spinorfield_cold";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_eoprec_spinorfield_cold";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "convert_from_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "saxpy";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "saxsbypz";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_zero_spinorfield";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "convert_to_kappa_format";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "convert_from_kappa_format";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "convert_to_kappa_format_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "convert_from_kappa_format_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "create_point_source";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "saxpy_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "saxsbypz_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "set_zero_spinorfield_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "create_point_source_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "scalar_product";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "scalar_product_reduction";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "global_squarenorm";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "global_squarenorm_reduction";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "scalar_product_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "global_squarenorm_eoprec";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "ratio";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "product";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "ps_correlator";
	Opencl::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	
}
#endif
