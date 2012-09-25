#include "opencl_module_correlator.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"
#include "meta/util.hpp"

using namespace std;

void Opencl_Module_Correlator::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Spinors::fill_collect_options(collect_options);
	//CP: give kappa and its negative value
	hmc_float kappa_tmp = get_parameters().get_kappa();
	*collect_options << " -DKAPPA=" << kappa_tmp;
	*collect_options << " -DMKAPPA=" << -kappa_tmp;

	if(get_parameters().get_use_pointsource() == true)
		*collect_options << " -DNUM_SOURCES=" << 12;
	else
		*collect_options << " -DNUM_SOURCES=" << get_parameters().get_num_sources();

	return;
}


void Opencl_Module_Correlator::fill_buffers()
{

	Opencl_Module_Spinors::fill_buffers();

	return;
}

void Opencl_Module_Correlator::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();
	basic_correlator_code = basic_fermion_code;
	logger.debug() << "Create correlator kernels...";

	if(get_parameters().get_use_pointsource() == true)
		create_point_source = createKernel("create_point_source") << basic_fermion_code << "spinorfield_point_source.cl";
	else
		create_stochastic_source = createKernel("create_stochastic_source") << basic_fermion_code << "spinorfield_stochastic_source.cl";

	switch (get_parameters().get_corr_dir()) {
		case 0 :
			correlator_ps = createKernel("correlator_ps_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_sc = createKernel("correlator_sc_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_vx = createKernel("correlator_vx_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_vy = createKernel("correlator_vy_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_vz = createKernel("correlator_vz_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_ax = createKernel("correlator_ax_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_ay = createKernel("correlator_ay_t") << basic_fermion_code << "fermionobservables.cl";
			correlator_az = createKernel("correlator_az_t") << basic_fermion_code << "fermionobservables.cl";
			break;
		case 3 :
			correlator_ps = createKernel("correlator_ps_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_sc = createKernel("correlator_sc_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_vx = createKernel("correlator_vx_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_vy = createKernel("correlator_vy_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_vz = createKernel("correlator_vz_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_ax = createKernel("correlator_ax_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_ay = createKernel("correlator_ay_z") << basic_fermion_code << "fermionobservables.cl";
			correlator_az = createKernel("correlator_az_z") << basic_fermion_code << "fermionobservables.cl";
			break;
		default:
			stringstream errmsg;
			errmsg << "Could not create correlator kernel as correlator direction " << get_parameters().get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}

	return;
}

void Opencl_Module_Correlator::clear_kernels()
{
	Opencl_Module_Spinors::clear_kernels();
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
	if(create_stochastic_source) {
		clerr = clReleaseKernel(create_stochastic_source);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	return;
}

void Opencl_Module_Correlator::clear_buffers()
{
	Opencl_Module_Spinors::clear_buffers();
	return;
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

void Opencl_Module_Correlator::create_point_source_device(cl_mem inout, int i, int spacepos, int timepos)
{
	set_zero_spinorfield_device(inout);
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_point_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_point_source, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 1, sizeof(int), &i);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 2, sizeof(int), &spacepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(create_point_source, 3, sizeof(int), &timepos);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( create_point_source, gs2, ls2);
}

void Opencl_Module_Correlator::create_stochastic_source_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(create_stochastic_source, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(create_stochastic_source, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	throw Opencl_Error(clerr, "stochastic source not yet implemented!!", __FILE__, __LINE__);
	get_device()->enqueue_kernel( create_stochastic_source, gs2, ls2);
}

void Opencl_Module_Correlator::correlator_device(const cl_kernel correlator_kernel, cl_mem in, cl_mem correlator)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(correlator_kernel, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(correlator_kernel, 0, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(correlator_kernel, 1, sizeof(cl_mem), &correlator);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(correlator_kernel , gs2, ls2);
}

#ifdef _PROFILING_
usetimer* Opencl_Module_Correlator::get_timer(const std::string& in) const
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Spinors::get_timer(in);
	if(noop != NULL) return noop;

	if (in == "create_point_source") {
		return &this->timer_create_point_source;
	}
	if (in == "create_stochastic_source") {
		return &this->timer_create_stochastic_source;
	}

	//CP: only one direction is calculated at a time, therefore one can use the same timer for both directions
	if (in == "correlator_ps_z") {
		return &this->timer_correlator_ps;
	}
	if (in == "correlator_sc_z") {
		return &this->timer_correlator_sc;
	}
	if (in == "correlator_vx_z") {
		return &this->timer_correlator_vx;
	}
	if (in == "correlator_vy_z") {
		return &this->timer_correlator_vy;
	}
	if (in == "correlator_vz_z") {
		return &this->timer_correlator_vz;
	}
	if (in == "correlator_ax_z") {
		return &this->timer_correlator_ax;
	}
	if (in == "correlator_ay_z") {
		return &this->timer_correlator_ay;
	}
	if (in == "correlator_az_z") {
		return &this->timer_correlator_az;
	}

	if (in == "correlator_ps_t") {
		return &this->timer_correlator_ps;
	}
	if (in == "correlator_sc_t") {
		return &this->timer_correlator_sc;
	}
	if (in == "correlator_vx_t") {
		return &this->timer_correlator_vx;
	}
	if (in == "correlator_vy_t") {
		return &this->timer_correlator_vy;
	}
	if (in == "correlator_vz_t") {
		return &this->timer_correlator_vz;
	}
	if (in == "correlator_ax_t") {
		return &this->timer_correlator_ax;
	}
	if (in == "correlator_ay_t") {
		return &this->timer_correlator_ay;
	}
	if (in == "correlator_az_t") {
		return &this->timer_correlator_az;
	}

	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}

#endif

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

#ifdef _PROFILING_

void Opencl_Module_Correlator::print_profiling(const std::string& filename, int number)
{
	Opencl_Module_Spinors::print_profiling(filename, number);
	if(create_point_source) {
		Opencl_Module::print_profiling(filename, create_point_source);
	}
	if(create_stochastic_source) {
		Opencl_Module::print_profiling(filename, create_stochastic_source);
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
#endif
