#include "opencl_module.hpp"

#include <fstream>

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"

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

const meta::Inputparameters& hardware::code::Opencl_Module::get_parameters() const noexcept
{
	return parameters;
}

hardware::Device * hardware::code::Opencl_Module::get_device() const noexcept
{
	return device;
}

TmpClKernel hardware::code::Opencl_Module::createKernel(const char * const kernel_name, std::string build_opts) const
{
	return device->create_kernel(kernel_name, build_opts);
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
	//CP: this is set manually to fit the longest kernel name
	out.width(32);
	out.precision(15);
	//to look like that
	/*
	logger.trace() << "*******************************************************************";
	logger.trace() << "Fermion\t"<< setfill(' ') << setw(16)<< "BW[GB/s]\t" << setfill(' ') << setw(18) << "Re/Wr[MByte]\t" << setfill(' ') << setw(6)  << "Calls\t" << setfill(' ') << setw(10)  << "Time[mus]";
	*/
	out << kernelName << "\t" << data.get_total_time() << "\t" << data.get_num_values() << "\t" << avg_time << "\t" << avg_time_site << "\t" << bandwidth << "\t" << flops << "\t" << (float) read_write_size / mega << "\t" << flop_size << std::endl;
	out.close();
	return;
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
	return;
}

void hardware::code::Opencl_Module::print_profiling(const std::string& filename, const cl_kernel& kernel) const
{
	// only print info if kernel has been initialized
	if(kernel) {
		const std::string kernel_name = get_kernel_name(kernel);
		const hardware::ProfilingData data = device->get_profiling_data(kernel);
		::print_profiling(filename, kernel_name, data, this->get_read_write_size(kernel_name), this->get_flop_size(kernel_name), meta::get_vol4d(get_parameters()));
	}
}

void hardware::code::Opencl_Module::print_profiling(const std::string& filename, int number) const
{
	logger.trace() << "Printing Profiling-information to file \"" << filename << "\"";
	print_profile_header(filename, number);
}
