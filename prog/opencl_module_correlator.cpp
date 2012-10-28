#include "opencl_module_correlator.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"
#include "meta/util.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device *, const meta::Inputparameters& params)
{
	std::ostringstream options;

	//CP: give kappa and its negative value
	hmc_float kappa_tmp = params.get_kappa();
	options << "-DKAPPA=" << kappa_tmp;
	options << " -DMKAPPA=" << -kappa_tmp;

	options << " -DNUM_SOURCES=" << params.get_num_sources();

	//CP: give content of sources as compile parameters
	options << " -DSOURCE_CONTENT=" << params.get_sourcecontent();

	return options.str();
}


void Opencl_Module_Correlator::fill_kernels()
{
	basic_correlator_code = basic_fermion_code << ClSourcePackage(collect_build_options(get_device(), get_parameters()));
	logger.debug() << "Create correlator kernels...";

	if(get_parameters().get_sourcetype() == meta::Inputparameters::point)
	  create_point_source = createKernel("create_point_source") << basic_correlator_code << prng_code << "spinorfield_point_source.cl";
	else if (get_parameters().get_sourcetype() == meta::Inputparameters::volume)
	  create_volume_source = createKernel("create_volume_source") << basic_correlator_code << prng_code << "spinorfield_volume_source.cl";
	else if (get_parameters().get_sourcetype() == meta::Inputparameters::timeslice)
	  create_timeslice_source = createKernel("create_timeslice_source") << basic_correlator_code << prng_code << "spinorfield_timeslice_source.cl";

	//CP: If a pointsource is chosen, the correlators have a particular simple form. 
	if(get_parameters().get_sourcetype() == meta::Inputparameters::point){
	  std::string filename_tmp =  "fermionobservables_correlators_point.cl";
	  switch (get_parameters().get_corr_dir()) {
	  case 0 :
	    correlator_ps = createKernel("correlator_ps_t") << basic_correlator_code << filename_tmp;
	    correlator_sc = createKernel("correlator_sc_t") << basic_correlator_code << filename_tmp;
	    correlator_vx = createKernel("correlator_vx_t") << basic_correlator_code << filename_tmp;
	    correlator_vy = createKernel("correlator_vy_t") << basic_correlator_code << filename_tmp;
	    correlator_vz = createKernel("correlator_vz_t") << basic_correlator_code << filename_tmp;
	    correlator_ax = createKernel("correlator_ax_t") << basic_correlator_code << filename_tmp;
	    correlator_ay = createKernel("correlator_ay_t") << basic_correlator_code << filename_tmp;
	    correlator_az = createKernel("correlator_az_t") << basic_correlator_code << filename_tmp;
	    break;
	  case 3 :
	    correlator_ps = createKernel("correlator_ps_z") << basic_correlator_code << filename_tmp;
	    correlator_sc = createKernel("correlator_sc_z") << basic_correlator_code << filename_tmp;
	    correlator_vx = createKernel("correlator_vx_z") << basic_correlator_code << filename_tmp;
	    correlator_vy = createKernel("correlator_vy_z") << basic_correlator_code << filename_tmp;
	    correlator_vz = createKernel("correlator_vz_z") << basic_correlator_code << filename_tmp;
	    correlator_ax = createKernel("correlator_ax_z") << basic_correlator_code << filename_tmp;
	    correlator_ay = createKernel("correlator_ay_z") << basic_correlator_code << filename_tmp;
	    correlator_az = createKernel("correlator_az_z") << basic_correlator_code << filename_tmp;
	    break;
	  default:
	    stringstream errmsg;
	    errmsg << "Could not create correlator kernel as correlator direction " << get_parameters().get_corr_dir() << " has not been implemented.";
	    throw Print_Error_Message(errmsg.str());
	  }
	} else{
	  std::string filename_tmp =  "fermionobservables_correlators_stochastic.cl";
	  switch (get_parameters().get_corr_dir()) {
	  case 0 :
	    correlator_ps = createKernel("correlator_ps_t") << basic_correlator_code << filename_tmp;
	    correlator_sc = createKernel("correlator_sc_t") << basic_correlator_code << filename_tmp;
	    correlator_vx = createKernel("correlator_vx_t") << basic_correlator_code << filename_tmp;
	    correlator_vy = createKernel("correlator_vy_t") << basic_correlator_code << filename_tmp;
	    correlator_vz = createKernel("correlator_vz_t") << basic_correlator_code << filename_tmp;
	    correlator_ax = createKernel("correlator_ax_t") << basic_correlator_code << filename_tmp;
	    correlator_ay = createKernel("correlator_ay_t") << basic_correlator_code << filename_tmp;
	    correlator_az = createKernel("correlator_az_t") << basic_correlator_code << filename_tmp;
	    break;
	  case 3 :
	    correlator_ps = createKernel("correlator_ps_z") << basic_correlator_code << filename_tmp;
	    correlator_sc = createKernel("correlator_sc_z") << basic_correlator_code << filename_tmp;
	    correlator_vx = createKernel("correlator_vx_z") << basic_correlator_code << filename_tmp;
	    correlator_vy = createKernel("correlator_vy_z") << basic_correlator_code << filename_tmp;
	    correlator_vz = createKernel("correlator_vz_z") << basic_correlator_code << filename_tmp;
	    correlator_ax = createKernel("correlator_ax_z") << basic_correlator_code << filename_tmp;
	    correlator_ay = createKernel("correlator_ay_z") << basic_correlator_code << filename_tmp;
	    correlator_az = createKernel("correlator_az_z") << basic_correlator_code << filename_tmp;
	    break;
	  default:
	    stringstream errmsg;
	    errmsg << "Could not create correlator kernel as correlator direction " << get_parameters().get_corr_dir() << " has not been implemented.";
	    throw Print_Error_Message(errmsg.str());
	  }
	}
}

void Opencl_Module_Correlator::clear_kernels()
{
	int clerr = CL_SUCCESS;
	clerr = clReleaseKernel(correlator_ps);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_sc);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_vx);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_vy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_vz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_ax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_ay);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(correlator_az);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(create_point_source) {
		clerr = clReleaseKernel(create_point_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(create_volume_source) {
		clerr = clReleaseKernel(create_volume_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(create_timeslice_source) {
		clerr = clReleaseKernel(create_timeslice_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
}

void Opencl_Module_Correlator::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module_Spinors::get_work_sizes(kernel, ls, gs, num_groups);

	//LZ: should be valid for all kernels for correlators, i.e. for names that look like correlator_??_?
	string kernelname = get_kernel_name(kernel);
	if( kernelname.find("correlator") == 0 ) {
		if(get_device()->get_device_type() == CL_DEVICE_TYPE_GPU) {
			*ls = get_parameters().get_nspace();
			*gs = *ls;
			*num_groups = 1;
		} else {
			*ls = 1;
			*gs = *ls;
			*num_groups = 1;
		}
	}

	return;
}

cl_kernel Opencl_Module_Correlator::get_correlator_kernel(string which)
{
	if( which.compare("ps") == 0 ) {
		return correlator_ps;
	}
	if( which.compare("sc") == 0 ) {
		return correlator_sc;
	}
	if( which.compare("vx") == 0 ) {
		return correlator_vx;
	}
	if( which.compare("vy") == 0 ) {
		return correlator_vy;
	}
	if( which.compare("vz") == 0 ) {
		return correlator_vz;
	}
	if( which.compare("ax") == 0 ) {
		return correlator_ax;
	}
	if( which.compare("ay") == 0 ) {
		return correlator_ay;
	}
	if( which.compare("az") == 0 ) {
		return correlator_az;
	}
	throw Print_Error_Message("get_correlator_kernel failed, no appropriate kernel found");
	return 0;
}

void Opencl_Module_Correlator::create_point_source_device(const hardware::buffers::Plain<spinor> * inout, int i, int spacepos, int timepos)
{
	set_zero_spinorfield_device(inout);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_point_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_point_source, gs2, ls2);

	if(logger.beDebug()) {
	  hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
	  hmc_float sqn;
	  this->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
	  sqn_tmp.dump(&sqn);
	  logger.debug() <<  "\t|source|^2:\t" << sqn;
	  if(sqn != sqn) {
	    throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}

}

void Opencl_Module_Correlator::create_volume_source_device(const hardware::buffers::Plain<spinor> * inout)
{
  set_zero_spinorfield_device(inout);

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_volume_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_volume_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_volume_source, 1, sizeof(cl_mem), get_prng_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_volume_source, gs2, ls2);

	if(logger.beDebug()) {
	  hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
	  hmc_float sqn;
	  this->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
	  sqn_tmp.dump(&sqn);
	  logger.debug() <<  "\t|source|^2:\t" << sqn;
	  if(sqn != sqn) {
	    throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}
}

void Opencl_Module_Correlator::create_timeslice_source_device(const hardware::buffers::Plain<spinor> * inout, const int timeslice)
{
  set_zero_spinorfield_device(inout);

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_timeslice_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_timeslice_source, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_timeslice_source, 1, sizeof(cl_mem), get_prng_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	int tmp = timeslice;
	clerr = clSetKernelArg(create_timeslice_source, 2, sizeof(int), &tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_timeslice_source, gs2, ls2);

	if(logger.beDebug()) {
	  hardware::buffers::Plain<hmc_float> sqn_tmp(1, get_device());
	  hmc_float sqn;
	  this->set_float_to_global_squarenorm_device(inout, &sqn_tmp);
	  sqn_tmp.dump(&sqn);
	  logger.debug() <<  "\t|source|^2:\t" << sqn;
	  if(sqn != sqn) {
	    throw Print_Error_Message("calculation of source gave nan! Aborting...", __FILE__, __LINE__);
	  }
	}
}

void Opencl_Module_Correlator::correlator_device(const cl_kernel correlator_kernel, const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<hmc_float> * correlator)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(correlator_kernel, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(correlator_kernel, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 1, sizeof(cl_mem), correlator->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(correlator_kernel , gs2, ls2);
}

size_t Opencl_Module_Correlator::get_read_write_size(const std::string& in) const
{
	size_t result = Opencl_Module_Spinors::get_read_write_size(in);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(parameters);
	//this returns the number of entries in an su3-matrix
	size_t R = meta::get_mat_size(parameters);
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (in == "create_point_source") {
		return 1000000000000000000000000;
	}
	if (in == "create_stochastic_source") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_ps_z" ) {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_sc_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vx_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vy_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_vz_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_ax_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_ay_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}
	if (in == "correlator_az_z") {
		//this kernel reads NUM_SOURCES spinors and writes NSPACE/NTIME real numbers
		int size_buffer = 0;
		int num_sources = get_parameters().get_num_sources();
		if(get_parameters().get_corr_dir() == 3) size_buffer = get_parameters().get_nspace();
		if(get_parameters().get_corr_dir() == 0) size_buffer = get_parameters().get_ntime();
		return num_sources * S * D * 12 * C + size_buffer * D;
	}

	return 0;
}

uint64_t Opencl_Module_Correlator::get_flop_size(const std::string& in) const
{
	uint64_t result = Opencl_Module_Spinors::get_flop_size(in);
	if (result != 0) return result;
	size_t S = meta::get_spinorfieldsize(get_parameters());
	size_t Seo = meta::get_eoprec_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	if (in == "create_point_source") {
		return 1000000000000000000000000;
	}
	if (in == "create_stochastic_source") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_ps_z" ) {
		return 1000000000000000000000000;
	}
	if (in == "correlator_sc_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_vx_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_vy_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_vz_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_ax_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_ay_z") {
		return 1000000000000000000000000;
	}
	if (in == "correlator_az_z") {
		return 1000000000000000000000000;
	}

	return 0;
}

void Opencl_Module_Correlator::print_profiling(const std::string& filename, int number)
{
	Opencl_Module_Spinors::print_profiling(filename, number);
	if(create_point_source) {
		Opencl_Module::print_profiling(filename, create_point_source);
	}
	if(create_volume_source) {
		Opencl_Module::print_profiling(filename, create_volume_source);
	}
	if(create_timeslice_source) {
		Opencl_Module::print_profiling(filename, create_timeslice_source);
	}
	Opencl_Module::print_profiling(filename, correlator_ps);
	Opencl_Module::print_profiling(filename, correlator_sc);
	Opencl_Module::print_profiling(filename, correlator_vx);
	Opencl_Module::print_profiling(filename, correlator_vy);
	Opencl_Module::print_profiling(filename, correlator_vz);
	Opencl_Module::print_profiling(filename, correlator_ax);
	Opencl_Module::print_profiling(filename, correlator_ay);
	Opencl_Module::print_profiling(filename, correlator_az);
}

Opencl_Module_Correlator::Opencl_Module_Correlator(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module_Spinors(params, device),
	  create_point_source(0), create_volume_source(0), create_timeslice_source(0),
	  correlator_ps(0), correlator_sc(0), correlator_vx(0), correlator_vy(0), correlator_vz(0), correlator_ax(0), correlator_ay(0), correlator_az(0)
{
	fill_kernels();
}

Opencl_Module_Correlator::~Opencl_Module_Correlator()
{
	clear_kernels();
}
