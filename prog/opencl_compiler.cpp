/** @file
 * Implementation of OpenCL compile process
 */

#include "opencl_compiler.hpp"

#include "logger.hpp"

#include <sstream>
#include <fstream>
#include <boost/regex.hpp>
#include <cstring>
#include <boost/algorithm/string.hpp>
#include "crypto/md5.h"

#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;
/// @todo quite some of this code could be simplified by moving from pure fstream to boost::filesystem equivalents

const std::string CACHE_DIR_NAME = "OpTiMaL/ocl_cache";

/**
 * Get the path on which to store a binary with the given md5
 *
 * This relies on boost to get the temporary directory, so boost methods
 * for specifying the temporary directory can be used.
 *
 * All filles will be stored in the subdirectory given by CACHE_DIR_NAME.
 */
static fs::path get_binary_file_path(std::string md5);

ClSourcePackage ClSourcePackage::operator <<(const std::string& file)
{
	auto tmp = this->files;
	tmp.push_back(file);
	return ClSourcePackage(tmp, this->options);
}

ClSourcePackage ClSourcePackage::operator <<(const ClSourcePackage& package)
{
	auto tmp = this->files;
	const auto other = package.files;
	tmp.insert(tmp.end(), other.begin(), other.end());
	return ClSourcePackage(tmp, options + ' ' + package.options);
}

ClSourcePackage ClSourcePackage::operator =(const ClSourcePackage& package)
{
	this->files = package.files;
	this->options = package.options;
	return *this;
}

const std::vector<std::string> ClSourcePackage::getFiles() const
{
	return files;
}

const std::string ClSourcePackage::getOptions() const
{
	return options;
}


TmpClKernel::operator cl_kernel() const
{
	cl_int clerr;

	logger.trace() << "Collecting sources to build the program for the " << kernel_name << " kernel";

	// write kernel files into sources
	// create array to point to contents of the different source files
	char ** sources = new char *[ files.size() ];
	size_t * source_sizes = new size_t[ files.size() ];

	std::string md5 = generateMD5();

	cl_program program = loadBinary(md5);

	if(!program) {
		logger.debug() << "Program not found in cache, building from source";
		std::string sourcecode;
		for(size_t n = 0; n < files.size(); n++) {
			std::stringstream tmp;
			tmp << SOURCEDIR << '/' << files[n];
			logger.debug() << "Read kernel source from file: " << tmp.str();

			std::fstream file;
			file.open(tmp.str().c_str());
			if( !file.is_open() ) throw File_Exception(tmp.str());

			file.seekg(0, std::ios::end);
			source_sizes[n] = file.tellg();
			file.seekg(0, std::ios::beg);

			sources[n] = new char[source_sizes[n]];

			file.read( sources[n], source_sizes[n] );

			file.close();
		}

		logger.trace() << "Creating program for the " << kernel_name << " kernel from collected sources";

		program = clCreateProgramWithSource(context, files.size() , (const char**) sources, source_sizes, &clerr);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateProgramWithSource", __FILE__, __LINE__);
	}

	logger.trace() << "Building kernel " << kernel_name << " using these options: \"" << build_options << "\"";

	clerr = clBuildProgram(program, 1, &device, build_options.c_str(), 0, 0);
	if(logger.beDebug()) {
		cl_int failed = clerr;
		if(clerr != CL_SUCCESS) {
			logger.error() << "... failed with error " << clerr << ", but look at BuildLog and abort then.";

			// dump program source
			size_t sourceSize;
			clerr = clGetProgramInfo(program, CL_PROGRAM_SOURCE, 0, NULL, &sourceSize);
			if(!clerr && sourceSize > 1) { // 0-terminated -> always at least one byte
				char* source = new char[sourceSize];
				clerr = clGetProgramInfo(program, CL_PROGRAM_SOURCE, sourceSize, source, &sourceSize);
				if(!clerr) {
					char const * const FILENAME = "broken_source.cl";
					std::ofstream srcFile(FILENAME);
					srcFile << source;
					srcFile.close();
					logger.debug() << "Dumped broken source to " << FILENAME;
				}
				delete[] source;
			}
		}

		logger.trace() << "Finished building program";

		// get build result
		size_t logSize;
		clerr = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
		if(!clerr && logSize > 1) { // 0-terminated -> always at least one byte
			logger.debug() << "Build Log:";
			char* log = new char[logSize];
			clerr = clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
			logger.debug() << log;
			delete [] log;
		}
		if(clerr) {
			throw Opencl_Error(clerr, "clGetProgramBuildInfo", __FILE__, __LINE__);
		}

		if(failed) {
			clerr = failed;
		}
	}
	if(clerr) {
		logger.fatal() << "... failed, aborting.";

		throw Opencl_Error(clerr, "clBuildProgram", __FILE__, __LINE__);
	}

	// store built binary to working directory
	dumpBinary(program, md5);

	// extract kernel
	cl_kernel kernel = clCreateKernel(program, kernel_name.c_str(), &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateKernel", __FILE__, __LINE__);

	if( logger.beDebug() ) {
		printResourceRequirements(kernel);
	}

	// make sure program get's cleaned up once kernel is released
	clReleaseProgram(program);

	// all done, hand back the kernel
	return kernel;
}

TmpClKernel TmpClKernel::operator <<(const std::string& file) const
{
	auto tmp = this->files;
	tmp.push_back(file);

	return TmpClKernel(kernel_name, build_options, context, device, tmp);
}

TmpClKernel TmpClKernel::operator <<(const ClSourcePackage& package) const
{
	auto tmp = this->files;
	const auto other = package.getFiles();
	tmp.insert(tmp.end(), other.begin(), other.end());

	return TmpClKernel(kernel_name, build_options + ' ' + package.getOptions(), context, device, tmp);
}

void TmpClKernel::printResourceRequirements(const cl_kernel kernel) const
{
	cl_int clerr;

	size_t nameSize;
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &nameSize );
	if( clerr == CL_SUCCESS ) {
		char* name = new char[nameSize];
		clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, nameSize, name, &nameSize );
		if( clerr == CL_SUCCESS )
			logger.trace() << "Kernel: " << name;
		delete[] name;
	}
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);

	// query the maximum work group size
	size_t work_group_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering maximum work group size is not supported on this device.";
		} else {
			throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
		}
	} else {
		logger.trace() << "  Maximum work group size: " << work_group_size;
	}

	// query the work group size specified at compile time (if any)
	size_t compile_work_group_size[3];
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3 * sizeof(size_t), compile_work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering compile time work group size is not supported on this device.";
		} else {
			throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
		}
	} else {
		if( compile_work_group_size[0] == 0 )
			logger.trace() << "  No work group size specified at compile time.";
		else
			logger.trace() << "  Compile time work group size: (" << compile_work_group_size[0] << ", " << compile_work_group_size[1] << ", " << compile_work_group_size[2] << ')';
	}

#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE // don't fail on OpenCL 1.0
	// query the preferred WORK_GROUP_SIZE_MULTIPLE (OpenCL 1.1 only)
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering work group size multiple is not supported on this device.";
		} else {
			throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
		}
	} else {
		logger.trace() << "  Preferred work group size multiple: " << work_group_size;
	}
#endif

	// query the local memory requirements
	cl_ulong local_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering local memory size is not supported on this device.";
		} else {
			throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
		}
	} else {
		logger.trace() << "  Local memory size (bytes): " << local_mem_size;
	}

#ifdef CL_KERNEL_PRIVATE_MEM_SIZE // don't fail on OpenCL 1.0
	// query the private memory required by the kernel (OpenCL 1.1 only)
	cl_ulong private_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(cl_ulong), &private_mem_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering private memory size is not supported on this device.";
		} else {
			throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
		}
	} else {
		logger.trace() << "  Private memory size (bytes): " << private_mem_size;
	}
#endif

	// the following only makes sense on AMD gpus ...
	cl_platform_id platform;
	clerr = clGetDeviceInfo( device, CL_DEVICE_PLATFORM, sizeof(cl_platform_id), &platform, NULL );
	if( clerr ) {
		logger.error() << "Failed to get the platform of the OpenCL device: ";
		return;
	}

	cl_device_type device_type;
	clerr = clGetDeviceInfo( device, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, NULL );
	if( clerr ) {
		logger.error() << "Failed to get the platform of the OpenCL device: ";
		return;
	}

	size_t platform_vendor_size;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 0, NULL, &platform_vendor_size);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		return;
	}
	char * platform_vendor = new char[platform_vendor_size];
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, platform_vendor_size, platform_vendor, NULL);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		return;
	}

	if( strcmp("Advanced Micro Devices, Inc.", platform_vendor) == 0
	    && device_type == CL_DEVICE_TYPE_GPU ) {

		// get device name
		size_t device_name_bytes;
		clerr = clGetDeviceInfo( device, CL_DEVICE_NAME, 0, NULL, &device_name_bytes );
		if( clerr ) {
			logger.error() << "Failed to get name of OpenCL device: ";
			return;
		}
		char * device_name = new char[device_name_bytes];
		clerr = clGetDeviceInfo( device, CL_DEVICE_NAME, device_name_bytes, device_name, NULL );
		if( clerr ) {
			logger.error() << "Failed to get name of OpenCL device: ";
			return;
		}

		logger.trace() << "Retrieving information for device " << device_name;

		size_t bytesInKernelName;
		clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName);
		if( clerr ) {
			logger.error() << "Failed to query kernel name: ";
			return;
		}
		char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
		clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL);
		if( clerr ) {
			logger.error() << "Failed to query kernel name: ";
			return;
		}

		logger.trace() << "Retrieving information for kernel " << kernelName;

		// retrieve some additinal info on the program
		std::stringstream tmp;
		tmp << kernelName << '_' << device_name << ".isa";
		std::string filename = tmp.str();

		logger.trace() << "Reading information from file " << filename;

		std::fstream isafile;
		isafile.open(filename.c_str());
		if(!isafile.is_open()) {
			logger.trace() << "Could not find ISA file under traditional name. Trying new name (Catalyst 12.4 and up).";
			std::stringstream tmp;
			tmp << "\\.\\/_temp_\\d+_" << device_name << '_' << kernelName << "\\.isa";
			boost::regex file_pattern(tmp.str(), boost::regex::icase);

			for(fs::directory_iterator i(fs::path(".")); i != fs::directory_iterator(); ++i) {
				std::string tmp = i->path().string();
				if(regex_match(tmp, file_pattern)) {
					isafile.open(tmp.c_str());
					break;
				}
			}

			if(!isafile.is_open()) {
				logger.error() << "Could not open ISA file. Aborting...";
				return;
			}
		}

		isafile.seekg(0, std::ios::end);
		size_t isasize = isafile.tellg();
		isafile.seekg(0, std::ios::beg);

		char * isabytes = new char[isasize];

		isafile.read( isabytes, isasize );

		isafile.close();

		std::string isa( isabytes );
		delete[] isabytes;

		boost::smatch what;

		// get scratch registers
		// starting with the Tahiti GPUs the format of the data changed,
		// we now also use the scratch register usage field as in indicator
		// for the file format
		boost::regex exScratch( "^MaxScratchRegsNeeded\\s*=\\s*(\\d*)$" );
		if( boost::regex_search( isa, what, exScratch ) ) {
			unsigned int scratch_regs = 0, gp_regs = 0, static_local_bytes = 0;

			logger.trace() << what[0];
			std::istringstream tmp( what[1] );
			tmp >> scratch_regs;

			// get GP registers
			boost::regex exGPR( "^SQ_PGM_RESOURCES:NUM_GPRS\\s*=\\s*(\\d*)$" );
			if( boost::regex_search( isa, what, exGPR ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> gp_regs;
			} else {
				logger.error() << "GPR usage section not found!";
			}

			// get GP registers
			boost::regex exStatic( "^SQ_LDS_ALLOC:SIZE\\s*=\\s*(0x\\d*)$" );
			if( boost::regex_search( isa, what, exStatic ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> std::hex >> static_local_bytes;
				static_local_bytes *= 4; // value in file is in units of floats
			} else {
				logger.error() << "Static local memory allocation section not found!";
			}

			logger.debug() << "Kernel: " << kernelName << " - " << gp_regs << " GPRs, " << scratch_regs << " scratch registers, "
			               << static_local_bytes << " bytes statically allocated local memory";
		} else {
			unsigned int scratch_regs = 0, sgp_regs = 0, vgp_regs = 0, static_local_bytes = 0;

			boost::regex exScratch( "^ScratchSize\\s*=\\s*(\\d*)\\s*;\\s*$" );
			if( boost::regex_search( isa, what, exScratch ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> scratch_regs;
			} else {
				logger.error() << "Scratch register usage section not found!";
			}

			// get GP registers
			boost::regex exSGPR( "^NumSgprs\\s*=\\s*(\\d*)\\s*;\\s*$" );
			if( boost::regex_search( isa, what, exSGPR ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> sgp_regs;
			} else {
				logger.error() << "sGPR usage section not found!";
			}

			// get GP registers
			boost::regex exVGPR( "^NumVgprs\\s*=\\s*(\\d*)\\s*;\\s*$" );
			if( boost::regex_search( isa, what, exVGPR ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> vgp_regs;
			} else {
				logger.error() << "vGPR usage section not found!";
			}

			// get GP registers
			boost::regex exStatic( "^COMPUTE_PGM_RSRC2:LDS_SIZE\s*=\s*(\d*)\s*$" );
			if( boost::regex_search( isa, what, exStatic ) ) {
				logger.trace() << what[0];
				std::istringstream tmp( what[1] );
				tmp >> static_local_bytes;
				static_local_bytes *= 4 * 64; // value in file is in units of floats
			} else {
				logger.trace() << "Static local memory allocation section not found. This is expected if no local memory is statically allocated.";
			}

			logger.debug() << "Kernel: " << kernelName << " - " << sgp_regs << " sGPRs, " << vgp_regs << " vGPRs, " << scratch_regs << " scratch registers, "
			               << static_local_bytes << " bytes statically allocated local memory";
		}
		delete[] device_name;
	} else {
		logger.trace() << "No AMD-GPU -> not scanning for kernel resource requirements";
	}

	delete[] platform_vendor;
}

std::string TmpClKernel::generateMD5() const
{
	// Points to consider for a unique binary identification
	//  * Device (Name + Driver Version): For simplicity assume we only compile for one device at a time (could also just add all devices)
	//  * Platform (Name + Version)
	//  * Compiler version (Should be given by platform and device driver version)
	//  * Sources names
	//  * Headers are sources, too
	//  * Build options

	/// @todo respect headers

	md5_t md5_state;
	md5_init(&md5_state);

	cl_int clerr;

	// add device info
	size_t device_name_bytes;
	clerr = clGetDeviceInfo( device, CL_DEVICE_NAME, 0, NULL, &device_name_bytes );
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL device: ";
		throw Opencl_Error(clerr);
	}
	char * device_name = new char[device_name_bytes];
	clerr = clGetDeviceInfo( device, CL_DEVICE_NAME, device_name_bytes, device_name, NULL );
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL device: ";
		throw Opencl_Error(clerr);
	}
	logger.trace() << "Adding " << device_name << " to MD5";
	md5_process(&md5_state, device_name, device_name_bytes);

	size_t driver_version_bytes;
	clerr = clGetDeviceInfo( device, CL_DRIVER_VERSION, 0, NULL, &driver_version_bytes );
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL device: ";
		throw Opencl_Error(clerr);
	}
	char * driver_version = new char[driver_version_bytes];
	clerr = clGetDeviceInfo( device, CL_DRIVER_VERSION, driver_version_bytes, driver_version, NULL );
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL device: ";
		throw Opencl_Error(clerr);
	}
	logger.trace() << "Adding " << driver_version << " to MD5";
	md5_process(&md5_state, driver_version, driver_version_bytes);

	// add platform information
	cl_platform_id platform;
	clerr = clGetDeviceInfo( device, CL_DEVICE_PLATFORM, sizeof(cl_platform_id), &platform, NULL );
	if( clerr ) {
		logger.error() << "Failed to get the platform of the OpenCL device: ";
		throw Opencl_Error(clerr);
	}

	size_t platform_name_bytes;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, 0, NULL, &platform_name_bytes);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		throw Opencl_Error(clerr);
	}
	char * platform_name = new char[platform_name_bytes];
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, platform_name_bytes, platform_name, NULL);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		throw Opencl_Error(clerr);
	}
	logger.trace() << "Adding " << platform_name << " to MD5";
	md5_process(&md5_state, platform_name, platform_name_bytes);

	size_t platform_version_bytes;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 0, NULL, &platform_version_bytes);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		throw Opencl_Error(clerr);
	}
	char * platform_version = new char[platform_version_bytes];
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, platform_version_bytes, platform_version, NULL);
	if( clerr ) {
		logger.error() << "Failed to get vendor of OpenCL platform: ";
		throw Opencl_Error(clerr);
	}
	logger.trace() << "Adding " << platform_version << " to MD5";
	md5_process(&md5_state, platform_version, platform_version_bytes);

	// add source information
	for(size_t n = 0; n < files.size(); n++) {
		md5_process(&md5_state, files[n].c_str(), files[n].length());
	}

	// build options
	logger.trace() << "Adding " << build_options << " to MD5";
	md5_process(&md5_state, build_options.c_str(), build_options.length());

	char sig[16];
	md5_finish(&md5_state, sig);

	char res[33];
	md5_sig_to_string(sig, res, 33);

	return std::string(res);
}

void TmpClKernel::dumpBinary(cl_program program, std::string md5) const
{

	cl_int clerr;

	// figure out which index is the one of the device we want to know
	cl_uint num_devices;
	clerr = clGetProgramInfo(program, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint), &num_devices, NULL);
	if( clerr ) {
		throw Opencl_Error(clerr);
	}

	cl_device_id * devices = new cl_device_id[num_devices];
	clerr = clGetProgramInfo(program, CL_PROGRAM_DEVICES, sizeof(cl_device_id) * num_devices, devices, NULL);
	if( clerr ) {
		throw Opencl_Error(clerr);
	}

	cl_uint device_index;
	for(device_index = 0; device_index < (num_devices - 1) && device != devices[device_index]; ++device_index) {
		// we only want to find the abortion criteria
	}

	size_t * program_binaries_bytes = new size_t[num_devices];
	clerr = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t) * num_devices, program_binaries_bytes, NULL);
	if( clerr ) {
		throw Opencl_Error(clerr);
	}
	size_t program_binary_bytes = program_binaries_bytes[device_index];
	if(program_binary_bytes == 0) {
		logger.warn() << "OpenCL implementation does not support exporting the compiled binary.";
		return;
	}

	unsigned char ** program_binaries = new unsigned char*[num_devices];
	unsigned char * program_binary = new unsigned char[program_binary_bytes];
	for(cl_uint i = 0; i < num_devices; ++i) {
		program_binaries[i] = (i == device_index) ? program_binary : NULL;
	}
	clerr = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(char*) * num_devices, program_binaries, NULL);
	if( clerr ) {
		throw Opencl_Error(clerr);
	}

	fs::path binaryfile = get_binary_file_path(md5);
	fs::create_directories(binaryfile.parent_path()); // make sure the directory exists, otherwise we cannot open the path for writing
	fs::ofstream outfile(binaryfile);
	outfile.write(reinterpret_cast<char*>(program_binary), program_binary_bytes);
	outfile.close();

	delete[] program_binary;
	delete[] program_binaries;
	delete[] program_binaries_bytes;
}

cl_program TmpClKernel::loadBinary(std::string md5) const
{
	/// @todo check whether any headers were modified
	fs::path binaryfile = get_binary_file_path(md5);
	logger.debug() << "Looking for cache file " << binaryfile;
	if(!fs::exists(binaryfile)) {
		return 0; // 0 is not a valid program object
	}

	std::time_t binary_date = fs::last_write_time(binaryfile);

	fs::path sourceDir(SOURCEDIR);
	for(size_t n = 0; n < files.size(); n++) {
		fs::path file = sourceDir / files[n];
		if(fs::last_write_time(file) > binary_date) {
			logger.debug() << "Sources have been modified since last program build. Recompilation is required.";
			return 0; // 0 is not a valid program object
		}
	}

	// read binary from file
	fs::ifstream inputfile(binaryfile);
	if(!inputfile.is_open()) {
		return 0; // 0 is not a valid program object
	}

	inputfile.seekg(0, std::ios::end);
	size_t binarysize = inputfile.tellg();
	inputfile.seekg(0, std::ios::beg);

	unsigned char * binary = new unsigned char[binarysize];

	inputfile.read( reinterpret_cast<char*>(binary), binarysize );

	inputfile.close();

	// create program from binary

	const unsigned char ** binaries = const_cast<const unsigned char**>(&binary); // we will only read from here on

	cl_int status, clerr;
	cl_program program = clCreateProgramWithBinary(context, 1, &device, &binarysize, binaries, &status, &clerr);
	if(clerr) {
		logger.error() << "Failed to reconstruct program from binary";
		throw Opencl_Error(clerr);
	}
	if(status) {
		logger.error() << "Failed to load program from binary";
		clReleaseProgram(program);
		program = 0;
	}

	delete[] binary;

	return program;
}


static fs::path get_binary_file_path(std::string md5)
{
	const fs::path cache_dir = fs::temp_directory_path() / CACHE_DIR_NAME;
	const std::string file_name = md5 + ".elf";
	return cache_dir / file_name;
}

TmpClKernel::TmpClKernel(const std::string kernel_name, const std::string build_options,
                         const cl_context context, cl_device_id device,
                         const std::vector<std::string> files)
	: kernel_name(kernel_name), build_options(boost::trim_copy(build_options)), context(context),
	  device(device), files(files) { }

ClSourcePackage::ClSourcePackage(const std::vector<std::string>& files, const std::string& options)
	: files(files), options(boost::trim_copy(options)) {}
