#include "opencl.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl::fill_collect_options(stringstream* collect_options)
{
	*collect_options << "-D_INKERNEL_ -DNSPACE=" << get_parameters()->get_ns() << " -DNTIME=" << get_parameters()->get_nt() << " -DVOLSPACE=" << get_parameters()->get_volspace();

	if(get_parameters()->get_use_rec12() == true)
		*collect_options << " -D_RECONSTRUCT_TWELVE_";
	if(get_parameters()->get_prec() == 64) {
		*collect_options << " -D_USEDOUBLEPREC_";
		if( device_double_extension.empty() ) {
			logger.warn() << "Warning: Undefined extension for use of double.";
		} else {
			*collect_options << " -D_DEVICE_DOUBLE_EXTENSION_" << device_double_extension << "_";
		}
	}
	if(get_parameters()->get_use_gpu() == true)
		*collect_options << " -D_USEGPU_";
	if(get_parameters()->get_use_chem_pot_re() == true) {
		*collect_options << " -D_CP_REAL_";
		*collect_options << " -DCPR=" << get_parameters()->get_chem_pot_re();
		*collect_options << " -DEXPCPR=" << exp(get_parameters()->get_chem_pot_re() );
		*collect_options << " -DMEXPCPR=" << exp(-1.*get_parameters()->get_chem_pot_re() );
	}
	if(get_parameters()->get_use_chem_pot_im() == true) {
		*collect_options << " -D_CP_IMAG_";
		*collect_options << " -DCPI=" << get_parameters()->get_chem_pot_im();
		*collect_options << " -DCOSCPI=" << cos( get_parameters()->get_chem_pot_im() );
		*collect_options << " -DSINCPI=" << sin( get_parameters()->get_chem_pot_im() );
	}
	if(get_parameters()->get_use_smearing() == true) {
		*collect_options << " -D_USE_SMEARING_";
		*collect_options << " -DRHO=" << get_parameters()->get_rho();
		*collect_options << " -DRHO_ITER=" << get_parameters()->get_rho_iter();
	}
	*collect_options << " -I" << SOURCEDIR;

	return;
}

cl_mem Opencl::get_clmem_gaugefield()
{
	return clmem_gaugefield;
}

cl_device_type Opencl::get_device_type()
{
	return device_type;
}

usetimer * Opencl::get_copy_on()
{
	return &copy_on;
}

usetimer * Opencl::get_copy_to()
{
	return &copy_to;
}

cl_mem Opencl::create_rw_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_READ_WRITE, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl::create_wo_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl::create_ro_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_READ_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl::create_uhp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_USE_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl::create_ahp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_ALLOC_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl::create_chp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

void Opencl::fill_buffers()
{
	logger.trace() << "Create buffer for gaugefield...";
	clmem_gaugefield = create_rw_buffer(NDIM * VOLSPACE * NTIME * sizeof(ocl_s_gaugefield));

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(sizeof(hmc_ocl_ran) * get_num_rndstates());

	logger.trace() << "Create buffer for gaugeobservables...";
	clmem_plaq = create_rw_buffer(sizeof(hmc_float));
	clmem_splaq = create_rw_buffer(sizeof(hmc_float));
	clmem_tplaq = create_rw_buffer(sizeof(hmc_float));
	clmem_polyakov = create_rw_buffer(sizeof(hmc_complex));

	// scratch buffers for gauge observable will be created on demand
	clmem_plaq_buf_glob = 0;
	clmem_tplaq_buf_glob = 0;
	clmem_splaq_buf_glob = 0;
	clmem_polyakov_buf_glob = 0;

	return;
}

void Opencl::fill_kernels()
{
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "opencl_geometry.cl" << "opencl_operations_complex.cl"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";


	logger.debug() << "Create gaugeobservables kernels...";
	plaquette = createKernel("plaquette") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	plaquette_reduction = createKernel("plaquette_reduction") << basic_opencl_code << "gaugeobservables_plaquette.cl";

	polyakov = createKernel("polyakov") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	polyakov_reduction = createKernel("polyakov_reduction") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	//only init if if wanted
	if(get_parameters()->get_use_smearing() == true) {
		stout_smear = createKernel("stout_smear") << basic_opencl_code << "gaugemomentum.cl" << "stout_smear.cl";
	}

}

void Opencl::init(cl_device_type wanted_device_type, inputparameters* params, int nstates)
{

	/** Number of threads to use for OpenCL kernels */
	if(params->get_use_gpu() == true) {
		numthreads = 128;
	} else {
		numthreads = 1;
	}

	init_basic(wanted_device_type, params, nstates);
	return;
}

void Opencl::init_basic(cl_device_type wanted_device_type, inputparameters* params, int nstates)
{
	init_rndarray(nstates);

	//variables, initializing, ...
	set_parameters(params);
	cl_int clerr = CL_SUCCESS;

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
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetPlatformIDs", __FILE__, __LINE__);


	//Cout Platforminfo
	char info[512];
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetPlatformInfo", __FILE__, __LINE__);
	logger.info() << "\tCL_PLATFORM_NAME:     " << info;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetPlatformInfo", __FILE__, __LINE__);
	logger.info() << "\tCL_PLATFORM_VENDOR:   " << info;
	clerr = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetPlatformInfo", __FILE__, __LINE__);
	logger.info() << "\tCL_PLATFORM_VERSION:  " << info;

	//Initializing devices
	cl_uint num_devices;
	clerr = clGetDeviceIDs(platform, wanted_device_type, 0, NULL, &num_devices);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceIDs", __FILE__, __LINE__);
	if(num_devices == 1) {
		logger.info() << "\t" << num_devices << " device of wanted type has been found.";
	} else {
		logger.info() << "\t" << num_devices << " devices of wanted type have been found. Choosing device number " << 0 << ".";
	}
	clerr = clGetDeviceIDs(platform, wanted_device_type, 1, &device, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceIDs", __FILE__, __LINE__);

	logger.info() << "\tDevice information: ";
	clerr = clGetDeviceInfo(device, CL_DEVICE_NAME, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\tCL_DEVICE_NAME:    " << info;
	clerr = clGetDeviceInfo(device, CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\tCL_DEVICE_VENDOR:  " << info;
	clerr = clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	if(device_type == CL_DEVICE_TYPE_CPU) logger.info() << "\t\tCL_DEVICE_TYPE:    CPU";
	if(device_type == CL_DEVICE_TYPE_GPU) logger.info() << "\t\tCL_DEVICE_TYPE:    GPU";
	if(device_type == CL_DEVICE_TYPE_ACCELERATOR) logger.info() << "\t\tCL_DEVICE_TYPE:    ACCELERATOR";
	if(device_type != CL_DEVICE_TYPE_CPU && device_type != CL_DEVICE_TYPE_GPU && device_type != CL_DEVICE_TYPE_ACCELERATOR)
		throw Print_Error_Message("Unexpected CL_DEVICE_TYPE...", __FILE__, __LINE__);
	clerr = clGetDeviceInfo(device, CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\tCL_DEVICE_VERSION: " << info;
	clerr = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\tCL_DEVICE_EXTENSIONS: " << info;

	if( strstr( info, "cl_amd_fp64" ) != NULL ) device_double_extension = "AMD";
	if( strstr( info, "cl_khr_fp64" ) != NULL ) device_double_extension = "KHR";


	// figure out the number of "cores"
	clerr = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);;
	logger.info() << "\t\tCL_DEVICE_MAX_COMPUTE_UNITS: " << max_compute_units;

	//Initilize context
	logger.trace() << "Create context...";
	context = clCreateContext(0, 1, &device, 0, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateContext", __FILE__, __LINE__);

	//Initilize queue
	logger.trace() << "Create command queue...";
#ifdef _PROFILING_
	queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &clerr);
#else
	queue = clCreateCommandQueue(context, device, 0, &clerr);
#endif
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateCommandQueue", __FILE__, __LINE__);

	//Create buffer
	this->fill_buffers();

	//Create kernels
	this->fill_kernels();

	//finish
	set_init_true();
	return;
}

void Opencl::finalize()
{

	cl_int clerr = CL_SUCCESS;


	if(get_init_status() == 1) {
		clerr = clFlush(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFlush", __FILE__, __LINE__);
		clerr = clFinish(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);
		this->clear_kernels();

		this->clear_buffers();

		clerr = clReleaseCommandQueue(queue);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseCommandQueue", __FILE__, __LINE__);
		clerr = clReleaseContext(context);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseContext", __FILE__, __LINE__);

		set_init_false();
	}
	return;
}

void Opencl::clear_kernels()
{
	logger.trace() << "Clearing kernels";

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(plaquette);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(polyakov);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(plaquette_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(polyakov_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters()->get_use_smearing() == true) {
		clerr = clReleaseKernel(stout_smear);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}

	return;
}

void Opencl::clear_buffers()
{
	logger.trace() << "Clearing memory objects";

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseMemObject(clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	clerr = clReleaseMemObject(clmem_plaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	clerr = clReleaseMemObject(clmem_tplaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	clerr = clReleaseMemObject(clmem_splaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	clerr = clReleaseMemObject(clmem_polyakov);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	if(clmem_plaq_buf_glob) {
		clerr = clReleaseMemObject(clmem_plaq_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	if(clmem_tplaq_buf_glob) {
		clerr = clReleaseMemObject(clmem_tplaq_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	if(clmem_splaq_buf_glob) {
		clerr = clReleaseMemObject(clmem_splaq_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
	if(clmem_polyakov_buf_glob) {
		clerr = clReleaseMemObject(clmem_polyakov_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}

	return;
}

void Opencl::copy_gaugefield_to_device(Matrixsu3* gaugefield)
{
	(*this->get_copy_to()).reset();
	const size_t gaugefield_size = NDIM * VOLSPACE * NTIME * sizeof(ocl_s_gaugefield);
	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(gaugefield_size);

	copy_to_ocl_format(host_gaugefield, gaugefield);

	cl_int clerr = clEnqueueWriteBuffer(queue, clmem_gaugefield, CL_TRUE, 0, gaugefield_size, host_gaugefield, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	free(host_gaugefield);

	(*this->get_copy_to()).add();
	return;
}


void Opencl::get_gaugefield_from_device(Matrixsu3* gaugefield)
{
	(*this->get_copy_to()).reset();
	const size_t gaugefield_size = NDIM * VOLSPACE * NTIME * sizeof(ocl_s_gaugefield);
	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(gaugefield_size);

	cl_int clerr = clEnqueueReadBuffer(queue, clmem_gaugefield, CL_TRUE, 0, gaugefield_size, host_gaugefield, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	copy_from_ocl_format(gaugefield, host_gaugefield);

	free(host_gaugefield);

	(*this->get_copy_to()).add();
	return;
}

void Opencl::copy_rndarray_to_device(hmc_ocl_ran* rndarray)
{
	(*this->get_copy_to()).reset();

	cl_int clerr = clEnqueueWriteBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_ocl_ran) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
	return;
}

void Opencl::copy_rndarray_from_device(hmc_ocl_ran* rndarray)
{
	(*this->get_copy_to()).reset();

	cl_int clerr = clEnqueueReadBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_ocl_ran) * get_num_rndstates(), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
	return;
}

void Opencl::copy_buffer_on_device(cl_mem in, cl_mem out, size_t size)
{
	(*this->get_copy_on()).reset();
	int clerr = clEnqueueCopyBuffer(queue, in, out, 0, 0, size , 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);

	(*this->get_copy_on()).add();
}

void Opencl::copy_buffer_to_device(void * source, cl_mem dest, size_t size)
{
	(*this->get_copy_to()).reset();

	int clerr = clEnqueueWriteBuffer(queue, dest, CL_TRUE, 0, size, source, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl::get_buffer_from_device(cl_mem source, void * dest, size_t size)
{
	(*this->get_copy_to()).reset();
	cl_int clerr = clEnqueueReadBuffer(queue, source, CL_TRUE, 0, size, dest, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size)
{
	///@todo make this properly handle multiple dimensions
	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = Opencl::get_numthreads(); /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sens on CPU

	// query the work group size specified at compile time (if any)
	size_t compile_work_group_size[3];
	cl_int clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3 * sizeof(size_t), compile_work_group_size, NULL );
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);

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

		delete [] kernelName;

	}

	cl_int clerr_enqueue = CL_SUCCESS;
#ifdef _PROFILING_
	cl_event event;
	clerr_enqueue = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, &event); //clerr error handling below
	if(clerr_enqueue == CL_SUCCESS) {

		cl_int done = clWaitForEvents(1, &event);
		if(done != CL_SUCCESS) throw Opencl_Error(clerr, "clWaitForEvents", __FILE__, __LINE__);

		//CP: Now I have to get the right timer, called timer_"kernelname"
		//First Method: Construct the explicit timername
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
//  char timerName[7] = ("timer_");
//  strcat(timerName, kernelName);
//  //Problem: How to call the memory object?
//  (this->timerName).add(get_kernel_exec_time(event));

		//Second Method: Nasty workaround
		//noop is used in case the kernel is not recognized
		usetimer *noop = NULL;
		noop = Opencl::get_timer(kernelName);
		if(noop == NULL)
			logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
		else
			(*get_timer(kernelName)).add(get_kernel_exec_time(event));

		delete [] kernelName;

	}
#else
	clerr_enqueue = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, NULL);
#endif
	if(clerr_enqueue != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting...";
		logger.fatal() << "Some more information:";
		logger.fatal() << "global_work_size: " << global_work_size;
		logger.fatal() << "local_work_size:  " << *local_work_size_p;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.fatal() << "Failed kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}

		throw Opencl_Error(clerr, "clEnqueueNDRangeKernel", __FILE__, __LINE__);

	}
}

void Opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size, const size_t local_work_size)
{
	cl_int clerr = CL_SUCCESS;
	cl_int clerr_enqueue = CL_SUCCESS;
#ifdef _PROFILING_
	cl_event event;
	clerr_enqueue = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, &event); //clerr evaluated below
	if(clerr_enqueue == CL_SUCCESS) {
		int done = clWaitForEvents(1, &event);
		if(done != CL_SUCCESS) throw Opencl_Error(clerr, "clWaitForEvents", __FILE__, __LINE__);

		//CP: Now I have to get the right timer, called timer_"kernelname"
		//First Method: Construct the explicit timername
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
		//  char timerName[7] = ("timer_");
//  strcat(timerName, kernelName);
//  //Problem: How to call the memory object?
//  (this->timerName).add(get_kernel_exec_time(event));

		//Second Method: Nasty workaround
		//noop is used in case the kernel is not recognized
		usetimer *noop = NULL;
//  noop = Opencl::get_timer(kernelName);
		noop = get_timer(kernelName);
		if(noop == NULL)
			logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
		else
			(*get_timer(kernelName)).add(get_kernel_exec_time(event));

		delete [] kernelName;
	}
#else
	clerr_enqueue = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
#endif

	if(clerr_enqueue != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting...";
		logger.fatal() << "Some more information:";
		logger.fatal() << "global_work_size: " << global_work_size;
		logger.fatal() << "local_work_size:  " << local_work_size;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.fatal() << "Failed kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}

		throw Opencl_Error(clerr, "clEnqueueNDRangeKernel", __FILE__, __LINE__);
	}
}
void Opencl::set_init_true()
{
	isinit = 1;
	return;
}
void Opencl::set_init_false()
{
	isinit = 0;
	return;
}
int Opencl::get_init_status()
{
	return isinit;
}

void Opencl::set_parameters (inputparameters * parameters_val)
{
	parameters = parameters_val;
	return;
}

inputparameters * Opencl::get_parameters ()
{
	return  parameters;
}

void Opencl::printResourceRequirements(const cl_kernel kernel)
{
	cl_int clerr = CL_SUCCESS;

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
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
	logger.trace() << "  Maximum work group size: " << work_group_size;

	// query the work group size specified at compile time (if any)
	size_t compile_work_group_size[3];
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_COMPILE_WORK_GROUP_SIZE, 3 * sizeof(size_t), compile_work_group_size, NULL );
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);

	if( compile_work_group_size[0] == 0 )
		logger.trace() << "  No work group size specified at compile time.";
	else
		logger.trace() << "  Compile time work group size: (" << compile_work_group_size[0] << ", " << compile_work_group_size[1] << ", " << compile_work_group_size[2] << ')';

#ifdef CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE // don't fail on OpenCL 1.0
	// query the preferred WORK_GROUP_SIZE_MULTIPLE (OpenCL 1.1 only)
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(size_t), &work_group_size, NULL );
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
	logger.trace() << "  Preferred work group size multiple: " << work_group_size;
#endif

	// query the local memory requirements
	cl_ulong local_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL );
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
	logger.trace() << "  Local memory size (bytes): " << local_mem_size;

#ifdef CL_KERNEL_PRIVATE_MEM_SIZE // don't fail on OpenCL 1.0
	// query the private memory required by the kernel (OpenCL 1.1 only)
	cl_ulong private_mem_size;
	clerr = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_PRIVATE_MEM_SIZE, sizeof(cl_ulong), &private_mem_size, NULL );
	if( clerr != CL_SUCCESS ) throw Opencl_Error(clerr, "clGetKernelWorkGroupInfo", __FILE__, __LINE__);
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

void Opencl::plaquette_device(cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes2(plaquette, this->get_device_type(), &ls, &gs, &num_groups);

	logger.debug() << "init scratch buffers if not already done";
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;

	if( clmem_plaq_buf_glob == 0 ) clmem_plaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_tplaq_buf_glob == 0 ) clmem_tplaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_splaq_buf_glob == 0 ) clmem_splaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_polyakov_buf_glob == 0 ) clmem_polyakov_buf_glob = create_rw_buffer(global_buf_size_complex);

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local plaquette calculation and first part of reduction
	int clerr = clSetKernelArg(plaquette, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(plaquette, gs, ls);

	// run second part of plaquette reduction

	this->get_work_sizes2(plaquette_reduction, this->get_device_type(), &ls, &gs, &num_groups);

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), &clmem_plaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), &clmem_tplaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), &clmem_splaq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	enqueueKernel(plaquette_reduction, gs, ls);

}

void Opencl::polyakov_device(cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes2(polyakov, this->get_device_type(), &ls, &gs, &num_groups);
	int buf_loc_size_complex = sizeof(hmc_complex) * ls;

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(polyakov, gs, ls);

	// second part of polyakov reduction

	this->get_work_sizes2(polyakov_reduction, this->get_device_type(), &ls, &gs, &num_groups);

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), &clmem_polyakov);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	///@todo improve
	ls = 1;
	gs = 1;
	enqueueKernel(polyakov_reduction, gs, ls);

}

void Opencl::gaugeobservables(cl_mem gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	//measure plaquette
	plaquette_device(gf);

	//read out values
	hmc_float plaq = 0.;
	hmc_float splaq = 0.;
	hmc_float tplaq = 0.;
	//NOTE: these are blocking calls!
	get_buffer_from_device(clmem_plaq, &plaq, sizeof(hmc_float));
	get_buffer_from_device(clmem_tplaq, &tplaq, sizeof(hmc_float));
	get_buffer_from_device(clmem_splaq, &splaq, sizeof(hmc_float));

	tplaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1));
	splaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	plaq  /= static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) * NC) / 2.;

	(*plaq_out) = plaq;
	(*splaq_out) = splaq;
	(*tplaq_out) = tplaq;

	//measure polyakovloop
	polyakov_device(gf);

	//read out values
	hmc_complex pol = hmc_complex_zero;
	//NOTE: this is a blocking call!
	get_buffer_from_device(clmem_polyakov, &pol, sizeof(hmc_complex));

	pol.re /= static_cast<hmc_float>(NC * VOLSPACE);
	pol.im /= static_cast<hmc_float>(NC * VOLSPACE);

	pol_out->re = pol.re;
	pol_out->im = pol.im;
}

TmpClKernel Opencl::createKernel(const char * const kernel_name)
{
	stringstream collect_options;
	this->fill_collect_options(&collect_options);
	return TmpClKernel(kernel_name, collect_options.str(), context, &device, 1);
}

void Opencl::stout_smear_device()
{

	return;
}

void Opencl::get_work_sizes2(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	//Construct explicit kernel name
	int clerr;
	size_t bytesInKernelName;
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);
	char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);

	/// @todo use kernelname
	size_t local_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		local_work_size = Opencl::get_numthreads(); /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sense on CPU

	size_t global_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * Opencl::get_numthreads() * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;

	//write out values
	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;

	delete [] kernelName;

	return;
}

void Opencl::get_work_sizes(size_t * ls, size_t * gs, cl_uint * num_groups, cl_device_type dev_type, string name)
{
	/// @todo use kernelname
	size_t local_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		local_work_size = Opencl::get_numthreads(); /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sense on CPU

	size_t global_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * Opencl::get_numthreads() * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;

	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;

	return;
}

#ifdef _PROFILING_
usetimer* Opencl::get_timer(char * in)
{
	if (strcmp(in, "polyakov_reduction") == 0) {
		return &(this->timer_polyakov_reduction);
	}
	if (strcmp(in, "polyakov") == 0) {
		return &(this->timer_polyakov);
	}
	if (strcmp(in, "plaquette_reduction") == 0) {
		return &(this->timer_plaquette_reduction);
	}
	if (strcmp(in, "plaquette") == 0) {
		return &(this->timer_plaquette);
	}
	if (strcmp(in, "stout_smear") == 0) {
		return &(this->timer_stout_smear);
	}
	//if the kernelname has not matched, return NULL
	else {
		return NULL;
	}
}

int Opencl::get_read_write_size(char * in, inputparameters * parameters)
{
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
		S = EOPREC_SPINORFIELDSIZE;
	else
		S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "polyakov") == 0) {
		return VOL4D * D * R + 1;
	}
	if (strcmp(in, "polyakov_reduction") == 0) {
		//this is not right, since one does not know bufelements now
		//return (Bufel + 1) *2
		return Opencl::get_numthreads();
	}
	if (strcmp(in, "plaquette") == 0) {
		return 48 * VOL4D * D * R + 1;
	}
	if (strcmp(in, "plaquette_reduction") == 0) {
		//this is not right, since one does not know bufelements now
		//return (Bufel + 1) *2
		return Opencl::get_numthreads();
	}
	if (strcmp(in, "stout_smear") == 0) {
		return 1000000000000000000000;
	}
	return 0;
}

void Opencl::print_profiling(std::string filename, char * kernelName, uint64_t time_total, int calls_total, int read_write_size)
{
	hmc_float bandwidth = 0.;
	uint64_t avg_time = 0.;
	uint64_t avg_time_site = 0.;
	//check if kernel has been called at all
	if(calls_total != 0 && time_total != 0) {
		avg_time = (uint64_t) ( ( (float) time_total ) / ((float) calls_total) );
		avg_time_site = (uint64_t) ( ( (float) time_total ) / ((float) (calls_total * VOL4D)) );
		//Bandwidth in GB/s: 1e-3 = 1e6 (museconds) * 1e-9 (GByte)
		bandwidth = (hmc_float) read_write_size / (hmc_float) time_total * (hmc_float) calls_total * 1e-3;
	}
	float mega = 1024 * 1024;
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	out.width(8);
	out.precision(15);
	//to look like that
	/*
	logger.trace() << "*******************************************************************";
	logger.trace() << "Fermion\t"<< setfill(' ') << setw(16)<< "BW[GB/s]\t" << setfill(' ') << setw(18) << "Re/Wr[MByte]\t" << setfill(' ') << setw(6)  << "Calls\t" << setfill(' ') << setw(10)  << "Time[mus]";
	*/
	out << kernelName << "\t" << time_total << "\t" << calls_total << "\t" << avg_time << "\t" << avg_time_site << "\t" << bandwidth << "\t" << (float) read_write_size / mega << std::endl;
	out.close();
	return;
}

void Opencl::print_profiling(std::string filename)
{
	logger.trace() << "Printing Profiling-information to file \"" << filename << "\"";
	char * kernelName;
	kernelName = "polyakov";
	print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "polyakov_reduction";
	print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "plaquette";
	print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "plaquette_reduction";
	print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
	kernelName = "stout_smear";
	print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters) );
}
#endif

void Opencl::init_rndarray(int nstates)
{
	num_rndstates = nstates;
}


int Opencl::get_numthreads()
{
	return numthreads;
}

int Opencl::get_num_rndstates()
{
	return num_rndstates;
}

