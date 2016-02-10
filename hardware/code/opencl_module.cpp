/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "opencl_module.hpp"

#include "../size_4.hpp"

using namespace std;

static void print_profile_header(const std::string& filename, int number);

/**
 * Print the profiling information of a specific kernel to a file.
 *
 * @param filename Name of file where data is appended.
 * @param kernelName Name of specific kernel.
 * @param time_total total execution time
 * @param calls_total total number of kernel calls
 * @param read_write_size number of bytes read and written by the kernel
 * @param flop_size amount of flops performed by the kernel
 */
static void print_profiling(const std::string& filename, const std::string& kernelName, const hardware::ProfilingData& data, size_t read_write_size, uint64_t flop_size, uint64_t sites);

static std::string collect_build_options(const hardware::Device * device, const hardware::code::OpenClKernelParametersInterface &kernelParameters)
{
	/**
	 * @Todo: Move all the explicit lattice size fcts. to own fcts.
	 * @Todo: Remove the "using"
	 */
	using namespace hardware::buffers;
	using namespace hardware::code;

	const size_4 local_size = device->getLocalLatticeExtents();
	const size_4 mem_size = device->getLocalLatticeMemoryExtents();

	std::ostringstream options;
	options.precision(16);
	
	options << "-I " << SOURCEDIR;
	options << " -D _INKERNEL_";
	options << " -D NSPACE=" << kernelParameters.getNs();

	options << " -D NTIME_GLOBAL=" << kernelParameters.getNt();
	options << " -D NTIME_LOCAL=" << local_size.t;
	options << " -D NTIME_MEM=" << mem_size.t;
	options << " -D NTIME_OFFSET=" << device->getGridPos().t * local_size.t;

	options << " -D VOLSPACE=" << kernelParameters.getSpatialLatticeVolume();

	options << " -D VOL4D_GLOBAL=" << kernelParameters.getLatticeVolume();
	options << " -D VOL4D_LOCAL=" << get_vol4d(local_size);
	options << " -D VOL4D_MEM=" << get_vol4d(mem_size);

	//this is needed for hmc_ocl_su3matrix
	options << " -D SU3SIZE=" << NC*NC << " -D STAPLEMATRIXSIZE=" << NC*NC;

	if(kernelParameters.getPrecision() == 64) {
		options << " -D _USEDOUBLEPREC_";
		// TODO renable support for older AMD GPUs
		//if( device_double_extension.empty() ) {
		//  logger.warn() << "Warning: Undefined extension for use of double.";
		//} else {
		//  options << " -D _DEVICE_DOUBLE_EXTENSION_" << device_double_extension << "_";
		//}
		options << " -D _DEVICE_DOUBLE_EXTENSION_KHR_";
	}
	if(device->get_device_type() == CL_DEVICE_TYPE_GPU )
		options << " -D _USEGPU_";
	if(kernelParameters.getUseChemPotRe() == true) {
		options << " -D _CP_REAL_";
		options << " -D CPR=" << kernelParameters.getChemPotRe();
		options << " -D EXPCPR=" << exp(kernelParameters.getChemPotRe() );
		options << " -D MEXPCPR=" << exp(-1.*kernelParameters.getChemPotRe() );
	}
	if(kernelParameters.getUseChemPotIm() == true) {
		options << " -D _CP_IMAG_";
		options << " -D CPI=" << kernelParameters.getChemPotIm();
		options << " -D COSCPI=" << cos( kernelParameters.getChemPotIm() );
		options << " -D SINCPI=" << sin( kernelParameters.getChemPotIm() );
	}
	if(kernelParameters.getUseSmearing() == true) {
		options << " -D _USE_SMEARING_";
		options << " -D RHO=" << kernelParameters.getRho();
		options << " -D RHO_ITER=" << kernelParameters.getRhoIter();
	}
	if(device->get_prefers_soa()) {
		options << " -D _USE_SOA_";
	}
	if(check_SU3_for_SOA(device)) {
		options << " -D GAUGEFIELD_STRIDE=" << get_SU3_buffer_stride(get_vol4d(mem_size) * NDIM, device);
	}
	if(check_Matrix3x3_for_SOA(device)) {
		options << " -D GAUGEFIELD_3X3_STRIDE=" << get_Matrix3x3_buffer_stride(get_vol4d(mem_size) * NDIM, device);
	}
	if(device->get_prefers_blocked_loops()) {
		options << " -D _USE_BLOCKED_LOOPS_";
	}
	if(kernelParameters.getUseRectangles() == true) {
		options <<  " -D _USE_RECT_" ;
		options <<  " -D C0=" << kernelParameters.getC0() << " -D C1=" << kernelParameters.getC1();
	}
	if(kernelParameters.getUseRec12() == true) {
		options <<  " -D _USE_REC12_" ;
	}
	if(kernelParameters.getUseEo()) {
		options << " -D EOPREC_SPINORFIELDSIZE_GLOBAL=" << kernelParameters.getEoprecSpinorFieldSize();
		options << " -D EOPREC_SPINORFIELDSIZE_LOCAL=" << get_vol4d(local_size)/2;
		options << " -D EOPREC_SPINORFIELDSIZE_MEM=" << get_vol4d(mem_size)/2;
	}
	//Always have non eo-prec options (for example in Wilson non eo kernels are always built)
	options << " -D SPINORFIELDSIZE_GLOBAL=" << kernelParameters.getSpinorFieldSize();
	options << " -D SPINORFIELDSIZE_LOCAL=" << get_vol4d(local_size);
	options << " -D SPINORFIELDSIZE_MEM=" << get_vol4d(mem_size);
	
	options << " -D GAUGEMOMENTASIZE_GLOBAL=" << kernelParameters.getLatticeVolume() * NDIM;
	options << " -D GAUGEMOMENTASIZE_LOCAL=" << get_vol4d(local_size) * NDIM;
	options << " -D GAUGEMOMENTASIZE_MEM=" << get_vol4d(mem_size) * NDIM;
	
	if(check_Gaugemomentum_for_SOA(device)) {
		options << " -D GAUGEMOMENTA_STRIDE=" << get_Gaugemomentum_buffer_stride(get_vol4d(mem_size)*NDIM, device);
	}
	if(check_Spinor_for_SOA(device)) {
		options << " -D EOPREC_SPINORFIELD_STRIDE=" << get_Spinor_buffer_stride(get_vol4d(mem_size)/2, device);
	}
	if(check_su3vec_for_SOA(device)) {
		options << " -D EOPREC_SU3VECFIELD_STRIDE=" << get_su3vec_buffer_stride(get_vol4d(mem_size)/2, device);
	}

	switch (kernelParameters.getFermact()) {
		case common::action::twistedmass :
			options << " -D _TWISTEDMASS_";
			break;
		case common::action::clover :
			options << " -D _CLOVER_";
			break;
		case common::action::rooted_stagg :
			options << " -D _RHMC_";
			options << " -D RA_MAX_ORDER=" << std::max(kernelParameters.getMetroApproxOrd(), kernelParameters.getMdApproxOrd());
			break;
		case common::action::wilson :
		case common::action::tlsym :
		case common::action::iwasaki :
		case common::action::dbw2 :
			// nothing to add
			break;
	}

	//CP: These are the BCs in spatial and temporal direction
	hmc_float tmp_spatial = (kernelParameters.getThetaFermionSpatial() * PI) / ( (hmc_float) kernelParameters.getNs());
	hmc_float tmp_temporal = (kernelParameters.getThetaFermionTemporal() * PI) / ( (hmc_float) kernelParameters.getNt());
	//BC: on the corners in each direction: exp(i theta) -> on each site exp(i theta*PI /LATEXTENSION) = cos(tmp2) + isin(tmp2)
	options << " -D SPATIAL_RE=" << cos(tmp_spatial);
	options << " -D MSPATIAL_RE=" << -cos(tmp_spatial);
	options << " -D SPATIAL_IM=" << sin(tmp_spatial);
	options << " -D MSPATIAL_IM=" << -sin(tmp_spatial);

	options << " -D TEMPORAL_RE=" << cos(tmp_temporal);
	options << " -D MTEMPORAL_RE=" << -cos(tmp_temporal);
	options << " -D TEMPORAL_IM=" << sin(tmp_temporal);
	options << " -D MTEMPORAL_IM=" << -sin(tmp_temporal);

	//This is mainly for molecular dynamics
	options <<  " -D BETA=" << kernelParameters.getBeta();
	
	//Options for correlators
	if(kernelParameters.getFermact() != common::action::rooted_stagg){
		hmc_float kappa_tmp = kernelParameters.getKappa();
		options << " -D KAPPA=" << kappa_tmp;
		options << " -D MKAPPA=" << -kappa_tmp;
	}
	options << " -D NUM_SOURCES=" << kernelParameters.getNumSources();
	//CP: give content of sources as compile parameters
	options << " -D SOURCE_CONTENT=" << kernelParameters.getSourceContent();
	
	//Options for heatbath
	if(kernelParameters.getUseAniso() == true) {
		options << " -D _ANISO_";
		options << " -D XI_0=" << kernelParameters.getXi0();
	}
	
	return options.str();
}

static std::vector<std::string> collect_build_files()
{
	std::vector<std::string> out;
	out.push_back("opencl_header.cl");
	out.push_back("globaldefs.h");
	out.push_back("types.h");

	return out;
}

hardware::code::Opencl_Module::Opencl_Module(const hardware::code::OpenClKernelParametersInterface &kernelParameters, const hardware::Device * deviceIn):
		kernelParameters(&kernelParameters),
		device(deviceIn),
		basic_sources ( ClSourcePackage(collect_build_files(),
				collect_build_options(device, kernelParameters) )  )
{}

hardware::code::Opencl_Module::~Opencl_Module(){}

const hardware::Device * hardware::code::Opencl_Module::get_device() const noexcept
{
	return device;
}

ClSourcePackage hardware::code::Opencl_Module::get_basic_sources() const noexcept
{
	return basic_sources;
}

TmpClKernel hardware::code::Opencl_Module::createKernel(const char * const kernel_name, std::string build_opts) const
{
	return device->createKernel(kernel_name, build_opts);
}

void hardware::code::Opencl_Module::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	//Query kernel name
	string kernelname = get_kernel_name(kernel);

	size_t local_work_size = device->get_preferred_local_thread_num();
	size_t global_work_size = device->get_preferred_global_thread_num();

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;

	//write out values
	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;
}

string hardware::code::Opencl_Module::get_kernel_name(const cl_kernel kernel) const
{
	int clerr;
	size_t bytesInKernelName;
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);
	char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);

	string kernel_name(kernelName);
	delete [] kernelName;

	return kernel_name;
}

static void print_profiling(const std::string& filename, const std::string& kernelName, const hardware::ProfilingData& data, size_t read_write_size, uint64_t flop_size, uint64_t sites)
{
	hmc_float bandwidth = 0.;
	hmc_float flops = 0.;
	uint64_t avg_time = 0.;
	uint64_t avg_time_site = 0.;
	//check if kernel has been called at all
	if(data.get_num_values()) {
		avg_time = (uint64_t) ( ( (float) data.get_total_time() ) / ((float) data.get_num_values() ) );
		avg_time_site = (uint64_t) ( ( (float) data.get_total_time() ) / ((float) (data.get_num_values() * sites)) );
		//Bandwidth in GB/s: 1e-3 = 1e6 (museconds) * 1e-9 (GByte)
		bandwidth = (hmc_float) read_write_size / (hmc_float) data.get_total_time() * (hmc_float) data.get_num_values() * 1e-3;
		flops = (hmc_float) flop_size / (hmc_float) data.get_total_time() * (hmc_float) data.get_num_values() * 1e-3;
	}
	float mega = 1024 * 1024;
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	out.width(32);
	out.precision(15);
	out << kernelName << "\t" << data.get_total_time() << "\t" << data.get_num_values() << "\t" << avg_time << "\t" << avg_time_site << "\t" << bandwidth << "\t" << flops << "\t" << (float) read_write_size / mega << "\t" << flop_size << std::endl;
	out.close();
}

static void print_profile_header(const std::string& filename, int number)
{
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	//CP: this is set manually to fit the longest kernel name
	out.width(32);
	out.precision(15);
	out << "#device " << number << "\tTime [mus]\tCalls\tAvg Time [mus]\tAvg Time/Site [mus]\tBW [GB/s]\tFLOPS [GFLOP/s]\tRe/Wr [MB]\tFLOP" << std::endl;
}

void hardware::code::Opencl_Module::print_profiling(const std::string& filename, const cl_kernel& kernel) const
{
	// only print info if kernel has been initialized
	if(kernel) {
		const std::string kernel_name = get_kernel_name(kernel);
		const hardware::ProfilingData data = device->getProfilingData(kernel);
		::print_profiling(filename, kernel_name, data, this->get_read_write_size(kernel_name), this->get_flop_size(kernel_name), kernelParameters->getLatticeVolume());
	}
}

void hardware::code::Opencl_Module::print_profiling(const std::string& filename, int number) const
{
	logger.trace() << "Printing Profiling-information to file \"" << filename << "\"";
	print_profile_header(filename, number);
}
