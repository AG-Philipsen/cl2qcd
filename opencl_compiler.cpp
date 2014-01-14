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

/** @file
 * Implementation of OpenCL compile process
 */

#include "opencl_compiler.hpp"

/**
 * Get the path on which to store a binary with the given md5
 *
 * This relies on boost to get the temporary directory, so boost methods
 * for specifying the temporary directory can be used.
 *
 * All filles will be stored in the subdirectory given by CACHE_DIR_NAME.
 */
static fs::path get_binary_file_path(std::string md5);
/**
 * Get the lock for the binary file of the given md5.
 *
 * The lock might be stored in an extra file next to the binary.
 */
static ip::file_lock get_lock_file(std::string md5);
/**
 * Get the absolute path to a sourcefile of the given name.
 */
static fs::path get_source_file_path(std::string filename);

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
	using namespace ip;

	cl_int clerr;

	logger.trace() << "Collecting sources to build the program for the " << kernel_name << " kernel";

	std::string md5 = generateMD5();

	cl_program program;

	// Us the lock to ensure that the binary is not read and written at the same time
	file_lock lock_file = get_lock_file(md5);
	{
		sharable_lock<file_lock> lock_shared(lock_file);
		program = loadBinary(md5);
	}
	if(!program) {
		// yes, at this exact point the file is not locked. someone might create a working binary while we are waiting
		// we need a write lock while we write to the file.
		// also scopy reading of source and building into this
		// once we got the lock first check whether maybe someone else build a working binary in the meantime
		scoped_lock<file_lock> lock_exclusive(lock_file);
		program = loadBinary(md5);
		if(program) {
			// yeah, easy way out
			buildProgram(program);
		} else {
			// build and write (the writing is why we need the exclusive lock
			program = loadSources();
			buildProgram(program);
			dumpBinary(program, md5);
		}
	} else {
		buildProgram(program);
	}

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
			boost::regex exStatic( "^COMPUTE_PGM_RSRC2:LDS_SIZE\\s*=\\s*(\\d*)\\s*$" );
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

	std::ostringstream id_string;

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
	id_string << std::string{device_name, device_name_bytes};

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
	id_string << std::string{driver_version, driver_version_bytes};

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
	id_string << std::string{platform_name, platform_name_bytes};

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
	id_string << std::string{platform_version, platform_version_bytes};

	// add source information
	for(auto filename: files) {
		std::string abs_filename = get_source_file_path(filename).string();
		id_string << abs_filename;
	}

	// build options
	id_string << build_options;

	return crypto::md5(id_string.str());
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

	for(size_t n = 0; n < files.size(); n++) {
		fs::path file = get_source_file_path(files[n]);
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

cl_program TmpClKernel::loadSources() const
{
	// write kernel files into sources
	// create array to point to contents of the different source files
	char ** sources = new char *[ files.size() ];
	size_t * source_sizes = new size_t[ files.size() ];

	logger.debug() << "Program not found in cache, building from source";
	std::string sourcecode;
	for(size_t n = 0; n < files.size(); n++) {
		fs::path filename = get_source_file_path(files[n]);
		logger.debug() << "Read kernel source from file: " << filename;

		fs::ifstream file(filename);
		if( !file.is_open() ) throw File_Exception(filename.string());

		file.seekg(0, std::ios::end);
		source_sizes[n] = file.tellg();
		file.seekg(0, std::ios::beg);

		sources[n] = new char[source_sizes[n]];

		file.read( sources[n], source_sizes[n] );

		file.close();
	}

	logger.trace() << "Creating program for the " << kernel_name << " kernel from collected sources";

	cl_int clerr;
	cl_program program = clCreateProgramWithSource(context, files.size() , (const char**) sources, source_sizes, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateProgramWithSource", __FILE__, __LINE__);

	return program;
}

static fs::path get_binary_file_path(std::string md5)
{
	static std::string user_name;
	if(user_name.empty()) {
		char* _user_name = getenv("USER");
		if(!_user_name) {
			throw Print_Error_Message("Failed to get user name", __FILE__, __LINE__);
		}
		user_name = _user_name;
	}
	const fs::path cache_dir = fs::temp_directory_path() / (user_name + '-' + CACHE_DIR_NAME);
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

void TmpClKernel::buildProgram(cl_program program) const
{
	logger.trace() << "Building kernel " << kernel_name << " using these options: \"" << build_options << "\"";

	cl_int clerr = clBuildProgram(program, 1, &device, build_options.c_str(), 0, 0);
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
}

static ip::file_lock get_lock_file(std::string md5)
{
	using namespace fs;

	path lock_file_path = get_binary_file_path(md5).replace_extension("lock");
	std::string lock_file_name = lock_file_path.string();
	logger.debug() << "Using lock file " << lock_file_name;
	if(!exists(lock_file_path)) {
		create_directories(lock_file_path.parent_path());
		ofstream dummy(lock_file_path);
	}
	return ip::file_lock(lock_file_name.c_str());
}

std::string getDirectoryBasedOnExtension(std::string filename)
{
	const fs::path filenamePath(filename);
	auto fileExtension = filenamePath.extension();
	if( fileExtension.string().compare(KERNEL_EXTENSION) == 0)
	  {
	    return KERNEL_DIRECTORY;
	  }
	else if( fileExtension.string().compare(HEADER_EXTENSION) == 0 || fileExtension.string().compare(HEADER_EXTENSION2) == 0 )
	  {
	    return HEADER_DIRECTORY;
	  }
	else
	  {
	    std::string errorMessage = "The file " + filename + " does not seem to have an extension known to the opencl compiler. Aborting...";
	    throw Print_Error_Message(errorMessage, __FILE__, __LINE__);
	  }
}

static fs::path get_source_file_path(std::string filename)
{
	const fs::path sourceDir(SOURCEDIR);
	std::string directoryName = getDirectoryBasedOnExtension(filename);
	return sourceDir / directoryName / filename;
}
