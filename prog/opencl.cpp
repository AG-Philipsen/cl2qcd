#include "opencl.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

ClSourcePackage basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "opencl_geometry.cl" << "opencl_operations_complex.cl"
                                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

hmc_error Opencl::fill_kernels_file ()
{
	//give a list of all kernel-files
	//!!CP: LZ should update this
	cl_kernels_file.push_back("opencl_header.cl");
	cl_kernels_file.push_back("opencl_geometry.cl");
	cl_kernels_file.push_back("random.cl");
	cl_kernels_file.push_back("opencl_operations_complex.cl");
	cl_kernels_file.push_back("operations_matrix_su3.cl");
	cl_kernels_file.push_back("operations_matrix.cl");
	cl_kernels_file.push_back("operations_gaugefield.cl");
	return HMC_SUCCESS;
}

hmc_error Opencl::fill_collect_options(stringstream* collect_options)
{
	*collect_options << "-D_INKERNEL_ -DNSPACE=" << NSPACE << " -DNTIME=" << NTIME << " -DVOLSPACE=" << VOLSPACE;

	//CP: these have to match those in the cmake file
#ifdef _RECONSTRUCT_TWELVE_
	*collect_options << " -D_RECONSTRUCT_TWELVE_";
#endif
#ifdef _USEDOUBLEPREC_
	*collect_options << " -D_USEDOUBLEPREC_";
#endif
#ifdef _USEGPU_
	*collect_options << " -D_USEGPU_";
#endif
#ifdef _CP_REAL_
	*collect_options << " -D_CP_REAL_";
#endif
#ifdef _CP_IMAG_
	*collect_options << " -D_CP_IMAG_";
#endif
#ifdef _USE_SMEARING_
	*collect_options << " -D_USE_SMEARING_";
#endif
	*collect_options << " -I" << SOURCEDIR;

	return HMC_SUCCESS;
}


hmc_error Opencl::fill_buffers()
{
	cl_int clerr;

	logger.trace() << "Create buffer for gaugefield...";
	clmem_gaugefield = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(s_gaugefield), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_rndarray), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	logger.trace() << "Create buffer for gaugeobservables...";
	clmem_plaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_splaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_tplaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clmem_polyakov = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	// scratch buffers for gauge observable will be created on demand
	clmem_plaq_buf_glob = 0;
	clmem_tplaq_buf_glob = 0;
	clmem_splaq_buf_glob = 0;
	clmem_polyakov_buf_glob = 0;

	return HMC_SUCCESS;
}

hmc_error Opencl::fill_kernels(cl_program clprogram)
{
	//CP: this should only be done when the heatbath wants to be used!!
	if(get_parameters()->get_perform_heatbath() == 1) {

		logger.debug() << "Create heatbath kernels...";
		heatbath_even = createKernel("heatbath_even") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
		if( logger.beDebug() )
			printResourceRequirements( heatbath_even );
		heatbath_odd = createKernel("heatbath_odd") << basic_opencl_code << "random.cl" << "update_heatbath.cl";
		if( logger.beDebug() )
			printResourceRequirements( heatbath_odd );

		logger.debug() << "Create overrelax kernels...";
		overrelax_even = createKernel("overrelax_even") << basic_opencl_code << "random.cl" << "overrelax.cl";
		if( logger.beDebug() )
			printResourceRequirements( overrelax_even );
		overrelax_odd = createKernel("overrelax_odd") << basic_opencl_code << "random.cl" << "overrelax.cl";
		if( logger.beDebug() )
			printResourceRequirements( overrelax_odd );

	}

	logger.debug() << "Create gaugeobservables kernels...";
	plaquette = createKernel("plaquette") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	if( logger.beDebug() )
		printResourceRequirements( plaquette );
	plaquette_reduction = createKernel("plaquette_reduction") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	if( logger.beDebug() )
		printResourceRequirements( plaquette_reduction );

	polyakov = createKernel("polyakov") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	if( logger.beDebug() )
		printResourceRequirements( polyakov );
	polyakov_reduction = createKernel("polyakov_reduction") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	if( logger.beDebug() )
		printResourceRequirements( polyakov_reduction );

	return HMC_SUCCESS;
}

hmc_error Opencl::init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* params)
{
	hmc_error err = init_basic(wanted_device_type, timer, params);
	return err;
}

hmc_error Opencl::init_basic(cl_device_type wanted_device_type, usetimer* timer, inputparameters* params)
{
	//variables, initializing, ...
	set_parameters(params);
	hmc_error err;
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	// in debug scenarios make the compiler dump the compile results
	if( logger.beDebug() ) {
		// the cast is safe here, we don't want to modify the value later
		putenv(const_cast<char*>("GPU_DUMP_DEVICE_KERNEL=3"));
	}

	//Initialize OpenCL,
	logger.trace() << "OpenCL being initialized...";

	cl_uint num_platforms;
	//LZ: for now, stick to one platform without any further checks...
	clerr = clGetPlatformIDs(1, &platform, &num_platforms);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clGetPlatformIDs failed...";
		exit(HMC_OCLERROR);
	}

	//Cout Platforminfo
	char info[512];
	if(clGetPlatformInfo(platform, CL_PLATFORM_NAME, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\tCL_PLATFORM_NAME:     " << info;
	if(clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\tCL_PLATFORM_VENDOR:   " << info;
	if(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\tCL_PLATFORM_VERSION:  " << info;

	//Initializing devices
	cl_uint num_devices;
	clerr = clGetDeviceIDs(platform, wanted_device_type, 0, NULL, &num_devices);
	if(num_devices == 1) {
		logger.info() << "\t" << num_devices << " device of wanted type has been found.";
	} else {
		logger.info() << "\t" << num_devices << " devices of wanted type have been found. Choosing device number " << 0 << ".";
	}
	clerr = clGetDeviceIDs(platform, wanted_device_type, 1, &device, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clGetDeviceIDs failed...";
		exit(HMC_OCLERROR);
	}

	logger.info() << "\tDevice information: ";
	if(clGetDeviceInfo(device, CL_DEVICE_NAME, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\t\tCL_DEVICE_NAME:    " << info;
	if(clGetDeviceInfo(device, CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\t\tCL_DEVICE_VENDOR:  " << info;
	if(clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(device_type == CL_DEVICE_TYPE_CPU) logger.info() << "\t\tCL_DEVICE_TYPE:    CPU";
	if(device_type == CL_DEVICE_TYPE_GPU) logger.info() << "\t\tCL_DEVICE_TYPE:    GPU";
	if(device_type == CL_DEVICE_TYPE_ACCELERATOR) logger.info() << "\t\tCL_DEVICE_TYPE:    ACCELERATOR";
	if(device_type != CL_DEVICE_TYPE_CPU && device_type != CL_DEVICE_TYPE_GPU && device_type != CL_DEVICE_TYPE_ACCELERATOR) {
		logger.fatal() << "unexpected CL_DEVICE_TYPE...";
		exit(HMC_OCLERROR);
	}
	if(clGetDeviceInfo(device, CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\t\tCL_DEVICE_VERSION: " << info;
	if(clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\t\tCL_DEVICE_EXTENSIONS: " << info;

	// figure out the number of "cores"
	if(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	logger.info() << "\t\tCL_DEVICE_MAX_COMPUTE_UNITS: " << max_compute_units;

	//Initilize context
	logger.trace() << "Create context...";
	context = clCreateContext(0, 1, &device, 0, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	//Initilize queue
	logger.trace() << "Create command queue...";
	//	queue = clCreateCommandQueue(context, device, 0, &clerr);
	queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	//collect all kernel files
	err = this->fill_kernels_file();
	if( err )
		exit( HMC_OCLERROR );

	//write kernel files into sources
	// create array to point to contents of the different source files
	char ** sources = new char *[ cl_kernels_file.size() ];
	size_t * source_sizes = new size_t[ cl_kernels_file.size() ];

	string sourcecode;
	for(size_t n = 0; n < cl_kernels_file.size(); n++) {
		stringstream tmp;
		tmp << SOURCEDIR << '/' << cl_kernels_file[n];
		logger.debug() << "Read kernel source from file: " << tmp.str();

		fstream kernelsfile;
		kernelsfile.open(tmp.str().c_str());
		if(!kernelsfile.is_open()) {
			logger.fatal() << "Could not open file. Aborting...";
			exit(HMC_FILEERROR);
		}

		kernelsfile.seekg(0, ios::end);
		source_sizes[n] = kernelsfile.tellg();
		kernelsfile.seekg(0, ios::beg);

		sources[n] = new char[source_sizes[n]];

		kernelsfile.read( sources[n], source_sizes[n] );

		kernelsfile.close();
	}

	//Create program from sources
	logger.trace() << "Create program...";
	cl_program clprogram = clCreateProgramWithSource(context, cl_kernels_file.size() , (const char**) sources, source_sizes, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	logger.debug() << "Build program...";

	//account for options
	stringstream collect_options;
	this->fill_collect_options(&collect_options);
	string buildoptions = collect_options.str();
	logger.debug() << "\tbuild options:" << "\t" << buildoptions;

	clerr = clBuildProgram(clprogram, 1, &device, buildoptions.c_str(), 0, 0);
	if(clerr != CL_SUCCESS) {
		logger.error() << "... failed with error " << clerr << ", but look at BuildLog and abort then.";
	}

	logger.trace() << "finished building program";

	size_t logSize;
	clerr |= clGetProgramBuildInfo(clprogram, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
	if(logSize > 1 && logger.beDebug()) { // 0-terminated -> always at least one byte
		logger.debug() << "Build Log:";
		char* log = new char[logSize];
		clerr |= clGetProgramBuildInfo(clprogram, device, CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
		logger.debug() << log;
		delete [] log;
	}
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";

		// dump program source
		size_t sourceSize;
		clerr = clGetProgramInfo(clprogram, CL_PROGRAM_SOURCE, 0, NULL, &sourceSize);
		if(!clerr && sourceSize > 1 && logger.beDebug()) { // 0-terminated -> always at least one byte
			char* source = new char[sourceSize];
			clerr = clGetProgramInfo(clprogram, CL_PROGRAM_SOURCE, sourceSize, source, &sourceSize);
			if(!clerr) {
				char const * const FILENAME = "broken_source.cl";
				ofstream srcFile(FILENAME);
				srcFile << source;
				srcFile.close();
				logger.debug() << "Dumped broken source to " << FILENAME;
			}
			delete[] source;
		}

		exit(HMC_OCLERROR);
	}

	//Create buffer
	err = this->fill_buffers();
	if( err )
		exit( HMC_OCLERROR );

	//Create kernels
	err = this->fill_kernels(clprogram);
	if( err )
		exit( HMC_OCLERROR );

	// release program
	clReleaseProgram(clprogram);

	//finish
	set_init_true();
	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::finalize()
{
	if(get_init_status() == 1) {
		if(clFlush(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clFinish(queue) != CL_SUCCESS) exit(HMC_OCLERROR);

		this->clear_kernels();

		this->clear_buffers();

		if(clReleaseCommandQueue(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseContext(context) != CL_SUCCESS) exit(HMC_OCLERROR);

		set_init_false();
	}
	return HMC_SUCCESS;
}

hmc_error Opencl::clear_kernels()
{
	logger.trace() << "Clearing kernels";

	if(get_parameters()->get_perform_heatbath() == 1) {
		if(clReleaseKernel(heatbath_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(heatbath_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(overrelax_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(overrelax_odd) != CL_SUCCESS) exit(HMC_OCLERROR);
	}

	if(clReleaseKernel(plaquette) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(plaquette_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(polyakov_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);

	return HMC_SUCCESS;
}

hmc_error Opencl::clear_buffers()
{
	logger.trace() << "Clearing memory objects";

	if(clReleaseMemObject(clmem_gaugefield) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_rndarray) != CL_SUCCESS) exit(HMC_OCLERROR);

	if(clReleaseMemObject(clmem_plaq) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_tplaq) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_splaq) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseMemObject(clmem_polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clmem_plaq_buf_glob) if(clReleaseMemObject(clmem_plaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clmem_tplaq_buf_glob) if(clReleaseMemObject(clmem_tplaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clmem_splaq_buf_glob) if(clReleaseMemObject(clmem_splaq_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clmem_polyakov_buf_glob) if(clReleaseMemObject(clmem_polyakov_buf_glob) != CL_SUCCESS) exit(HMC_OCLERROR);

	return HMC_SUCCESS;
}

hmc_error Opencl::copy_gaugefield_to_device(s_gaugefield* gaugefield, usetimer* timer)
{
//   cout<<"Copy gaugefield to device..."<<endl;
	timer->reset();
	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(sizeof(s_gaugefield));

	copy_to_ocl_format(host_gaugefield, gaugefield);

	int clerr = clEnqueueWriteBuffer(queue, clmem_gaugefield, CL_TRUE, 0, sizeof(s_gaugefield), host_gaugefield, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "...copy gaugefield failed, aborting.";
		exit(HMC_OCLERROR);
	}

	free(host_gaugefield);

	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::copy_rndarray_to_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Copy randomarray to device..."<<endl;
	timer->reset();

	int clerr = clEnqueueWriteBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::get_gaugefield_from_device(s_gaugefield* gaugefield, usetimer* timer)
{
//   cout<<"Get gaugefield from device..."<<endl;
	timer->reset();
	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(sizeof(s_gaugefield));

	int clerr = clEnqueueReadBuffer(queue, clmem_gaugefield, CL_TRUE, 0, sizeof(s_gaugefield), host_gaugefield, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		logger.fatal() << "errorcode :" << clerr;
		exit(HMC_OCLERROR);
	}

	copy_from_ocl_format(gaugefield, host_gaugefield);

	free(host_gaugefield);

	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Get randomarray from device..."<<endl;
	timer->reset();

	int clerr = clEnqueueReadBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::run_heatbath(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;
	timer->reset();

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_even failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_even failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_even, global_work_size);
	}

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 at heatbath_odd failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg2 at heatbath_odd failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;

}

hmc_error Opencl::run_overrelax(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;

	timer->reset();

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
	else
		global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg4 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_even, global_work_size);
	}

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg7 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "clSetKernelArg8 failed, aborting...";
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;
}

hmc_error Opencl::gaugeobservables(hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, usetimer* timer1, usetimer * timer2)
{
	cl_int clerr = CL_SUCCESS;

	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU ) {
		// reductions are broken for local_work_size > 64
		local_work_size = 64;//NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	} else {
		local_work_size = 1; // nothing else makes sens on CPU
	}

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups;

	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;

	if( clmem_plaq_buf_glob == 0 ) {
		clmem_plaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating clmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_tplaq_buf_glob == 0 ) {
		clmem_tplaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating tclmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_splaq_buf_glob == 0 ) {
		clmem_splaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating sclmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_polyakov_buf_glob == 0 ) {
		clmem_polyakov_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_complex, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating clmem_polyakov_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}


	//measure plaquette
	timer1->reset();

	hmc_float plaq;
	hmc_float splaq;
	hmc_float tplaq;
	int buf_loc_size_float = sizeof(hmc_float) * local_work_size;
	int buf_loc_size_complex = sizeof(hmc_complex) * local_work_size;

	plaq = 0.;
	splaq = 0.;
	tplaq = 0.;

	// run local plaquette calculation and first part of reduction

	clerr = clSetKernelArg(plaquette, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg4 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette, global_work_size, local_work_size);

	// run second part of plaquette reduction

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), &clmem_plaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), &clmem_tplaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg4 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), &clmem_splaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_plaq, CL_FALSE, 0, sizeof(hmc_float), &plaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_tplaq, CL_FALSE, 0, sizeof(hmc_float), &tplaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_splaq, CL_FALSE, 0, sizeof(hmc_float), &splaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	// wait for results to have been read back
	clFinish(queue);

	//two plaquette-measurements per thread -> add. factor of 1/2
	tplaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1));
	splaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	plaq  /= static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) * NC) / 2.;

	(*plaq_out) = plaq;
	(*splaq_out) = splaq;
	(*tplaq_out) = tplaq;

	timer1->add();

	//measure polyakovloop
	timer2->reset();
	hmc_complex pol;
	pol = hmc_complex_zero;

	// local polyakov compuation and first part of reduction

	clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov, global_work_size, local_work_size);

	// second part of polyakov reduction

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), &clmem_polyakov);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_polyakov, CL_FALSE, 0, sizeof(hmc_complex), &pol, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	// wait for result to have been read back
	clFinish(queue);

	pol.re /= static_cast<hmc_float>(NC * VOLSPACE);
	pol.im /= static_cast<hmc_float>(NC * VOLSPACE);

	pol_out->re = pol.re;
	pol_out->im = pol.im;

	timer2->add();

	return HMC_SUCCESS;
}

void Opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size)
{
	///@todo make this properly handle multiple dimensions
	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sens on CPU

	// query the work group size specified at compile time (if any)
	size_t compile_work_group_size[3];
	cl_int clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3 * sizeof(size_t), compile_work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	const size_t * const local_work_size_p = (compile_work_group_size[0] == 0) ? &local_work_size : &compile_work_group_size[0];

	// make sure global_work_size is divisible by global_work_size
	if( global_work_size % *local_work_size_p ) {
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

		logger.fatal() << "Kernel " << kernelName << " can only be run with a global work size which is a multiple of " << *local_work_size_p << ". The requested size was " << global_work_size << '.';
	}

	clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting..." << clerr << " - " << global_work_size << " - " << *local_work_size_p;

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
		logger.fatal() << "Failed kernel: " << kernelName;

		exit(HMC_OCLERROR);
	}
}

void Opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size, const size_t local_work_size)
{
	cl_int clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting..." << clerr << " - " << global_work_size << " - " << local_work_size;

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
		logger.fatal() << "Failed kernel: " << kernelName;

		exit(HMC_OCLERROR);
	}
}
hmc_error Opencl::set_init_true()
{
	isinit = 1;
	return HMC_SUCCESS;
}
hmc_error Opencl::set_init_false()
{
	isinit = 0;
	return HMC_SUCCESS;
}
int Opencl::get_init_status()
{
	return isinit;
}

hmc_error Opencl::set_parameters (inputparameters * parameters_val)
{
	parameters = parameters_val;
	return HMC_SUCCESS;
}

inputparameters * Opencl::get_parameters ()
{
	return  parameters;
}

void Opencl::printResourceRequirements(const cl_kernel kernel)
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
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}

	// query the maximum work group size
	size_t work_group_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "  Maximum work group size: " << work_group_size;

	// query the work group size specified at compile time (if any)
	size_t compile_work_group_size[3];
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3 * sizeof(size_t), compile_work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	if( compile_work_group_size[0] == 0 )
		logger.trace() << "  No work group size specified at compile time.";
	else
		logger.trace() << "  Compile time work group size: (" << compile_work_group_size[0] << ", " << compile_work_group_size[1] << ", " << compile_work_group_size[2] << ')';

#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE // don't fail on OpenCL 1.0
	// query the preferred WORK_GROUP_SIZE_MULTIPLE (OpenCL 1.1 only)
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &work_group_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "  Preferred work group size multiple: " << work_group_size;
#endif

	// query the local memory requirements
	cl_ulong local_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "  Local memory size (bytes): " << local_mem_size;

#ifdef CL_KERNEL_PRIVATE_MEM_SIZE // don't fail on OpenCL 1.0
	// query the private memory required by the kernel (OpenCL 1.1 only)
	cl_ulong private_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(cl_ulong), &private_mem_size, NULL );
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "Querying kernel properties failed, aborting...";
		exit(HMC_OCLERROR);
	}
	logger.trace() << "  Private memory size (bytes): " << private_mem_size;
#endif

	// the following only makes sense on AMD gpus ...

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

//CP: this is the overloaded function to calc gaugeobservables with a specific gaugefield
hmc_error Opencl::gaugeobservables(cl_mem gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, usetimer* timer1, usetimer * timer2)
{
	cl_int clerr = CL_SUCCESS;

	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU ) {
		// reductions are broken for local_work_size > 64
		local_work_size = 64;//NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
	} else {
		local_work_size = 1; // nothing else makes sens on CPU
	}

	size_t global_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * NUMTHREADS * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups;

	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;

	if( clmem_plaq_buf_glob == 0 ) {
		clmem_plaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating clmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_tplaq_buf_glob == 0 ) {
		clmem_tplaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating tclmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_splaq_buf_glob == 0 ) {
		clmem_splaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating sclmem_plaq_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_polyakov_buf_glob == 0 ) {
		clmem_polyakov_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_complex, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			logger.fatal() << "creating clmem_polyakov_buf_glob failed, aborting...";
			exit(HMC_OCLERROR);
		}
	}


	//measure plaquette
	timer1->reset();

	hmc_float plaq;
	hmc_float splaq;
	hmc_float tplaq;
	int buf_loc_size_float = sizeof(hmc_float) * local_work_size;
	int buf_loc_size_complex = sizeof(hmc_complex) * local_work_size;

	plaq = 0.;
	splaq = 0.;
	tplaq = 0.;

	// run local plaquette calculation and first part of reduction

	clerr = clSetKernelArg(plaquette, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg4 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette, global_work_size, local_work_size);

	// run second part of plaquette reduction

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), &clmem_plaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg3 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), &clmem_tplaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg4 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), &clmem_splaq);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg5 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg6 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_plaq, CL_FALSE, 0, sizeof(hmc_float), &plaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_tplaq, CL_FALSE, 0, sizeof(hmc_float), &tplaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_splaq, CL_FALSE, 0, sizeof(hmc_float), &splaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	// wait for results to have been read back
	clFinish(queue);

	//two plaquette-measurements per thread -> add. factor of 1/2
	tplaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1));
	splaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	plaq  /= static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) * NC) / 2.;

	(*plaq_out) = plaq;
	(*splaq_out) = splaq;
	(*tplaq_out) = tplaq;

	timer1->add();

	//measure polyakovloop
	timer2->reset();
	hmc_complex pol;
	pol = hmc_complex_zero;

	// local polyakov compuation and first part of reduction

	clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov, global_work_size, local_work_size);

	// second part of polyakov reduction

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg0 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), &clmem_polyakov);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg1 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clSetKernelArg2 failed, aborting...";
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_polyakov, CL_FALSE, 0, sizeof(hmc_complex), &pol, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	// wait for result to have been read back
	clFinish(queue);

	pol.re /= static_cast<hmc_float>(NC * VOLSPACE);
	pol.im /= static_cast<hmc_float>(NC * VOLSPACE);

	pol_out->re = pol.re;
	pol_out->im = pol.im;

	timer2->add();

	return HMC_SUCCESS;
}

TmpClKernel Opencl::createKernel(const char * const kernel_name)
{
	stringstream collect_options;
	this->fill_collect_options(&collect_options);
	return TmpClKernel(kernel_name, collect_options.str(), context, &device, 1);
}
