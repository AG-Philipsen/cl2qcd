/** @file
 * Implementation of OpenCL compile process
 */

#include "opencl_compiler.hpp"

#include "logger.hpp"
#include "hmcerrs.h"

#include <sstream>
#include <fstream>

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
			logger.fatal() << "Could not open file. Aborting...";
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

	logger.trace() << "Building kernel " << kernel_name;

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

