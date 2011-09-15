#include "opencl_module.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

void Opencl_Module::init(cl_command_queue queue, cl_mem* clmem_gaugefield, inputparameters* params, int maxcomp, string double_ext)
{

	set_queue(queue);
	set_gaugefield(clmem_gaugefield);
	set_parameters(params);

	set_device_double_extension(double_ext);
	set_max_compute_units(maxcomp);

	cl_uint clerr = clGetCommandQueueInfo(get_queue(), CL_QUEUE_CONTEXT, sizeof(cl_context), &ocl_context, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetCommandQueueInfo", __FILE__, __LINE__);


	clerr = clGetCommandQueueInfo(get_queue(), CL_QUEUE_DEVICE, sizeof(cl_device_id), &device, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetCommandQueueInfo", __FILE__, __LINE__);

	clerr = clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);

	clerr = clGetDeviceInfo(device, CL_DEVICE_PLATFORM, sizeof(cl_platform_id), &platform, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);

	switch ( device_type ) {
		case CL_DEVICE_TYPE_GPU :
			numthreads = 128;
			break;
		case CL_DEVICE_TYPE_CPU :
			numthreads = 1;
			break;
		default :
			throw Print_Error_Message("Could not retrive proper CL_DEVICE_TYPE...", __FILE__, __LINE__);
	}

	this->fill_buffers();
	this->fill_kernels();
}


void Opencl_Module::finalize()
{

	this->clear_buffers();
	this->clear_kernels();
	return;
}


cl_context Opencl_Module::get_context()
{
	return ocl_context;
}


void Opencl_Module::set_queue(cl_command_queue queue)
{
	ocl_queue = queue;
	return;
}

cl_command_queue Opencl_Module::get_queue()
{
	return ocl_queue;
}

void Opencl_Module::set_gaugefield(cl_mem* clmem_gaugefield)
{
	ocl_gaugefield = clmem_gaugefield;
	return;
}

cl_mem* Opencl_Module::get_gaugefield()
{
	return ocl_gaugefield;
}

void Opencl_Module::set_parameters(inputparameters* params)
{
	parameters = params;
	return;
}

inputparameters* Opencl_Module::get_parameters()
{
	return parameters;
}

cl_device_type Opencl_Module::get_device_type()
{
	return device_type;
}

cl_device_id Opencl_Module::get_device()
{
	return device;
}

cl_platform_id Opencl_Module::get_platform()
{
	return platform;
}


void Opencl_Module::fill_collect_options(stringstream* collect_options)
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
	if( device_type == CL_DEVICE_TYPE_GPU )
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


cl_mem Opencl_Module::create_rw_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_READ_WRITE, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl_Module::create_wo_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_WRITE_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl_Module::create_ro_buffer(size_t size)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_READ_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl_Module::create_uhp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_USE_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl_Module::create_ahp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_ALLOC_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

cl_mem Opencl_Module::create_chp_buffer(size_t size, void *host_pointer)
{
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(get_context(), CL_MEM_COPY_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
	return tmp;
}

void Opencl_Module::fill_buffers()
{

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

void Opencl_Module::fill_kernels()
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

void Opencl_Module::clear_kernels()
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

void Opencl_Module::clear_buffers()
{
	logger.trace() << "Clearing memory objects";

	cl_int clerr = CL_SUCCESS;

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


void Opencl_Module::copy_buffer_on_device(cl_mem in, cl_mem out, size_t size)
{
	(*this->get_copy_on()).reset();

	int clerr = clEnqueueCopyBuffer(get_queue(), in, out, 0, 0, size , 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);

	(*this->get_copy_on()).add();
}

void Opencl_Module::copy_buffer_to_device(void * source, cl_mem dest, size_t size)
{
	(*this->get_copy_to()).reset();

	int clerr = clEnqueueWriteBuffer(get_queue(), dest, CL_TRUE, 0, size, source, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl_Module::get_buffer_from_device(cl_mem source, void * dest, size_t size)
{
	(*this->get_copy_to()).reset();
	cl_int clerr = clEnqueueReadBuffer(get_queue(), source, CL_TRUE, 0, size, dest, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl_Module::enqueueKernel(const cl_kernel kernel, const size_t global_work_size)
{
	///@todo make this properly handle multiple dimensions
	// decide on work-sizes
	size_t local_work_size;
	if( device_type == CL_DEVICE_TYPE_GPU )
		local_work_size = Opencl_Module::get_numthreads(); /// @todo have local work size depend on kernel properties (and device? autotune?)
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
	clerr_enqueue = clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, &event); //clerr error handling below
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
		noop = Opencl_Module::get_timer(kernelName);
		if(noop == NULL)
			logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
		else
			(*get_timer(kernelName)).add(get_kernel_exec_time(event));

		delete [] kernelName;

	}
#else
	clerr_enqueue = clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, NULL);
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

void Opencl_Module::enqueueKernel(const cl_kernel kernel, const size_t global_work_size, const size_t local_work_size)
{
	cl_int clerr = CL_SUCCESS;
	cl_int clerr_enqueue = CL_SUCCESS;
#ifdef _PROFILING_
	cl_event event;
	clerr_enqueue = clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, &event); //clerr evaluated below
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
//  noop = Opencl_Module::get_timer(kernelName);
		noop = get_timer(kernelName);
		if(noop == NULL)
			logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
		else
			(*get_timer(kernelName)).add(get_kernel_exec_time(event));

		delete [] kernelName;
	}
#else
	clerr_enqueue = clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
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

void Opencl_Module::printResourceRequirements(const cl_kernel kernel)
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
	clerr = clGetPlatformInfo(get_platform(), CL_PLATFORM_NAME, 0, NULL, &platform_name_size);
	if( clerr ) {
		logger.error() << "Failed to get name of OpenCL platform: ";
		return;
	}
	char * platform_name = new char[platform_name_size];
	clerr = clGetPlatformInfo(get_platform(), CL_PLATFORM_NAME, platform_name_size, platform_name, NULL);
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

void Opencl_Module::plaquette_device(cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(plaquette, this->get_device_type(), &ls, &gs, &num_groups);

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

	this->get_work_sizes(plaquette_reduction, this->get_device_type(), &ls, &gs, &num_groups);

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

void Opencl_Module::polyakov_device(cl_mem gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(polyakov, this->get_device_type(), &ls, &gs, &num_groups);
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

	this->get_work_sizes(polyakov_reduction, this->get_device_type(), &ls, &gs, &num_groups);

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

void Opencl_Module::gaugeobservables(cl_mem gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
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

TmpClKernel Opencl_Module::createKernel(const char * const kernel_name)
{
	stringstream collect_options;
	this->fill_collect_options(&collect_options);
	return TmpClKernel(kernel_name, collect_options.str(), get_context(), &device, 1);
}

void Opencl_Module::stout_smear_device()
{

	throw Print_Error_Message("Not implemented yet.");

	return;
}

void Opencl_Module::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	/// @todo use kernelname
	size_t local_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		local_work_size = Opencl_Module::get_numthreads(); /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sense on CPU

	size_t global_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * Opencl_Module::get_numthreads() * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;

	//write out values
	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;

	return;
}

string Opencl_Module::get_kernel_name(const cl_kernel kernel)
{
	int clerr;
	size_t bytesInKernelName;
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);
	char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);

	string kernel_name = kernelName;
	delete [] kernelName;

	return kernel_name;
}

usetimer * Opencl_Module::get_copy_on()
{
	return &copy_on;
}

usetimer * Opencl_Module::get_copy_to()
{
	return &copy_to;
}

#ifdef _PROFILING_
usetimer* Opencl_Module::get_timer(char * in)
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

int Opencl_Module::get_read_write_size(char * in, inputparameters * parameters)
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
		return Opencl_Module::get_numthreads();
	}
	if (strcmp(in, "plaquette") == 0) {
		return 48 * VOL4D * D * R + 1;
	}
	if (strcmp(in, "plaquette_reduction") == 0) {
		//this is not right, since one does not know bufelements now
		//return (Bufel + 1) *2
		return Opencl_Module::get_numthreads();
	}
	if (strcmp(in, "stout_smear") == 0) {
		return 1000000000000000000000;
	}
	return 0;
}

void Opencl_Module::print_profiling(std::string filename, const char * kernelName, uint64_t time_total, int calls_total, int read_write_size)
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

void Opencl_Module::print_profiling(std::string filename)
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


int Opencl_Module::get_numthreads()
{
	return numthreads;
}

void Opencl_Module::set_max_compute_units(int maxcomp)
{
	max_compute_units = maxcomp;
	return;
}

int Opencl_Module::get_max_compute_units()
{
	return max_compute_units;
}

void Opencl_Module::set_device_double_extension(string double_ext)
{
	if(double_ext.empty()) {
		device_double_extension.clear();
	} else {
		device_double_extension.assign(double_ext);
	}
	return;
}

string Opencl_Module::get_device_double_extension()
{
	return device_double_extension;
}

void Opencl_Module::print_copy_times(uint64_t totaltime)
{
	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time

	uint64_t copy2_time = (this->copy_on).getTime();
	uint64_t copy1_time = (this->copy_to).getTime();

	int copy1_steps =  (this->copy_to).getNumMeas();
	int copy2_steps =  (this->copy_on).getNumMeas();

	uint64_t copy1_avgtime = divide(copy1_time, copy1_steps);
	uint64_t copy2_avgtime = divide(copy2_time, copy2_steps);

	logger.trace() << "## *******************************************************************";
	logger.trace() << "## Copy-Times\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg" << '\t' << setw(5) << "perc";
	logger.trace() << "## CpyTo:\t" << setfill(' ') << setw(12) << copy1_time << '\t' << setw(12) << copy1_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(copy1_time, totaltime);
	logger.trace() << "## CpyOn:\t" << setfill(' ') << setw(12) << copy2_time << '\t' << setw(12) << copy2_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(copy2_time, totaltime);
	logger.trace() << "## *******************************************************************";

	logger.trace() << "## No output of times to file implemented yet...";
	/** @todo output to file is not implemented */
	//See older files for example code
	return;
}
