/** @file
 * Implementation of OpenCL compile process
 */

#include "opencl_compiler.hpp"

#include "logger.hpp"
#include "hmcerrs.h"

#include <sstream>
#include <fstream>
#include <boost/regex.hpp>
#include <cstring>

ClSourcePackage ClSourcePackage::operator <<(const char *file)
{
	std::vector<const char*> tmp = this->files;
	tmp.push_back(file);
	return ClSourcePackage(tmp);
}

ClSourcePackage ClSourcePackage::operator <<(const ClSourcePackage& package)
{
	std::vector<const char*> tmp = this->files;
	const std::vector<const char*> other = package.files;
	tmp.insert(tmp.end(), other.begin(), other.end());
	return ClSourcePackage(tmp);
}

ClSourcePackage ClSourcePackage::operator =(const ClSourcePackage& package)
{
	this->files = package.files;
	return *this;
}

const std::vector<const char *> ClSourcePackage::getFiles() const
{
	return files;
}


TmpClKernel::operator cl_kernel() const
{
	cl_int clerr;

	logger.trace() << "Collecting sources to build the program for the " << kernel_name << " kernel";

	//write kernel files into sources
	// create array to point to contents of the different source files
	char ** sources = new char *[ files.size() ];
	size_t * source_sizes = new size_t[ files.size() ];

	std::string sourcecode;
	for(size_t n = 0; n < files.size(); n++) {
		std::stringstream tmp;
		tmp << SOURCEDIR << '/' << files[n];
		logger.debug() << "Read kernel source from file: " << tmp.str();

		std::fstream file;
		file.open(tmp.str().c_str());
		if(!file.is_open()) {
			logger.fatal() << "Could not open file " << tmp.str() << ". Aborting...";
			exit(HMC_FILEERROR);
		}

		file.seekg(0, std::ios::end);
		source_sizes[n] = file.tellg();
		file.seekg(0, std::ios::beg);

		sources[n] = new char[source_sizes[n]];

		file.read( sources[n], source_sizes[n] );

		file.close();
	}

	logger.trace() << "Creating program for the " << kernel_name << " kernel from collected sources";

	cl_program program = clCreateProgramWithSource(context, files.size() , (const char**) sources, source_sizes, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	logger.trace() << "Building kernel " << kernel_name << " using these options: " << build_options;

	clerr = clBuildProgram(program, num_devices, devices, build_options.c_str(), 0, 0);
	if(clerr != CL_SUCCESS && logger.beDebug()) {
		logger.error() << "... failed with error " << clerr << ", but look at BuildLog and abort then.";
	}

	logger.trace() << "Finished building program";

	// get build result for each device
	for(size_t i = 0; i < num_devices; ++i) {
		size_t logSize;
		clerr |= clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
		if(logSize > 1 && logger.beDebug()) { // 0-terminated -> always at least one byte
			logger.debug() << "Build Log:";
			char* log = new char[logSize];
			clerr |= clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
			logger.debug() << log;
			delete [] log;
		}
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "... failed, aborting.";

			// dump program source
			size_t sourceSize;
			clerr = clGetProgramInfo(program, CL_PROGRAM_SOURCE, 0, NULL, &sourceSize);
			if(!clerr && sourceSize > 1 && logger.beDebug()) { // 0-terminated -> always at least one byte
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

			exit(HMC_OCLERROR);
		}
	}

	// extract kernel
	cl_kernel kernel = clCreateKernel(program, kernel_name, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Failed to create kernel from program.";
		exit(HMC_OCLERROR);
	}

	if( logger.beDebug() ) {
		for(size_t i = 0; i < num_devices; ++i)
			printResourceRequirements(kernel, devices[i]);
	}

	// make sure program get's cleaned up once kernel is released
	clReleaseProgram(program);

	// all done, hand back the kernel
	return kernel;
}

TmpClKernel TmpClKernel::operator <<(const char *file) const
{
	std::vector<const char*> tmp = this->files;
	tmp.push_back(file);

	return TmpClKernel(kernel_name, build_options, context, devices, num_devices, tmp);
}

TmpClKernel TmpClKernel::operator <<(const ClSourcePackage& package) const
{
	std::vector<const char*> tmp = this->files;
	const std::vector<const char*> other = package.getFiles();
	tmp.insert(tmp.end(), other.begin(), other.end());

	return TmpClKernel(kernel_name, build_options, context, devices, num_devices, tmp);
}

void TmpClKernel::printResourceRequirements(const cl_kernel kernel, const cl_device_id device) const
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
	if( clerr != CL_SUCCESS ) {
		logger.fatal() << "Querying kernel properties failed: " << clerr;
		exit(HMC_OCLERROR);
	}

	// query the maximum work group size
	size_t work_group_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		if( clerr == CL_INVALID_VALUE ) {
			logger.warn() << "Quering maximum work group size is not supported on this device.";
		} else {
			logger.fatal() << "Querying kernel properties failed: " << clerr;
			exit(HMC_OCLERROR);
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
			logger.fatal() << "Querying kernel properties failed: " << clerr;
			exit(HMC_OCLERROR);
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
			logger.fatal() << "Querying kernel properties failed: " << clerr;
			exit(HMC_OCLERROR);
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
			logger.fatal() << "Querying kernel properties failed: " << clerr;
			exit(HMC_OCLERROR);
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
			logger.fatal() << "Querying kernel properties failed: " << clerr;
			exit(HMC_OCLERROR);
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

	size_t platform_name_size;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, 0, NULL, &platform_name_size);
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL platform: ";
		return;
	}
	char * platform_name = new char[platform_name_size];
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, platform_name_size, platform_name, NULL);
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL platform: ";
		return;
	}

	if( strcmp("AMD Accelerated Parallel Processing", platform_name) == 0
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
			logger.error() << "Could not open ISA file. Aborting...";
			return;
		}

		isafile.seekg(0, std::ios::end);
		size_t isasize = isafile.tellg();
		isafile.seekg(0, std::ios::beg);

		char * isabytes = new char[isasize];

		isafile.read( isabytes, isasize );

		isafile.close();

		std::string isa( isabytes );
		delete[] isabytes;

		unsigned int scratch_regs, gp_regs, static_local_bytes;

		boost::smatch what;

		// get scratch registers
		boost::regex exScratch( "^MaxScratchRegsNeeded\\s*=\\s*(\\d*)$" );
		if( boost::regex_search( isa, what, exScratch ) ) {
			logger.trace() << what[0];
			std::istringstream tmp( what[1] );
			tmp >> scratch_regs;
		} else {
			logger.error() << "Scratch register usage section not found!";
		}

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

		delete[] device_name;
	} else {
		logger.trace() << "No AMD-GPU -> not scanning for kernel resource requirements";
	}

	delete[] platform_name;
}
