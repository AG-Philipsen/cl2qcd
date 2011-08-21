#include "opencl.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;

hmc_error Opencl::fill_collect_options(stringstream* collect_options)
{
	*collect_options << "-D_INKERNEL_ -DNSPACE=" << get_parameters()->get_nspace() << " -DNTIME=" << get_parameters()->get_ntime() << " -DVOLSPACE=" << get_parameters()->get_volspace();
	
	if(get_parameters()->get_use_rec12() == true)
		*collect_options << " -D_RECONSTRUCT_TWELVE_";
	if(get_parameters()->get_prec() == 64){
		*collect_options << " -D_USEDOUBLEPREC_";
		if( device_double_extension.empty() ) {
			logger.warn() << "Warning: Undefined extension for use of double.";
		} else {
			*collect_options << " -D_DEVICE_DOUBLE_EXTENSION_"<<device_double_extension<<"_";
		}
	}
	if(get_parameters()->get_use_gpu() == true)
		*collect_options << " -D_USEGPU_";
	if(get_parameters()->get_use_chem_pot_re() == true)
		*collect_options << " -D_CP_REAL_";
	if(get_parameters()->get_use_chem_pot_im() == true)
		*collect_options << " -D_CP_IMAG_";
	if(get_parameters()->get_use_smearing() == true){
		*collect_options << " -D_USE_SMEARING_";
		*collect_options << " -DRHO="<< get_parameters()->get_rho();
		*collect_options << " -DRHO_ITER="<< get_parameters()->get_rho_iter();
	}
	*collect_options << " -I" << SOURCEDIR;

	return HMC_SUCCESS;
}

cl_mem Opencl::get_clmem_gaugefield(){
	return clmem_gaugefield;
}

cl_device_type Opencl::get_device_type(){
	return device_type;
}

cl_mem Opencl::create_rw_buffer(size_t size){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_READ_WRITE, size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating read-write buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

cl_mem Opencl::create_wo_buffer(size_t size){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating write-only buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

cl_mem Opencl::create_ro_buffer(size_t size){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_READ_ONLY, size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating read-only buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

cl_mem Opencl::create_uhp_buffer(size_t size, void *host_pointer){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_USE_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating use-host-pointer buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

cl_mem Opencl::create_ahp_buffer(size_t size, void *host_pointer){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_ALLOC_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating alloc-host-pointer buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

cl_mem Opencl::create_chp_buffer(size_t size, void *host_pointer){
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(context, CL_MEM_COPY_HOST_PTR, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... creating copy-host-pointer buffer failed, aborting.";
		exit(HMC_OCLERROR);
	}
	return tmp;
}

hmc_error Opencl::fill_buffers()
{
	logger.trace() << "Create buffer for gaugefield...";
	clmem_gaugefield = create_rw_buffer(sizeof(s_gaugefield));

	logger.trace() << "Create buffer for random numbers...";
	clmem_rndarray = create_rw_buffer(sizeof(hmc_rndarray));

	logger.trace() << "Create buffer for gaugeobservables...";
	clmem_plaq = create_rw_buffer(sizeof(hmc_float) * global_work_size);
	clmem_splaq = create_rw_buffer(sizeof(hmc_float) * global_work_size);
	clmem_tplaq = create_rw_buffer(sizeof(hmc_float) * global_work_size);
	clmem_polyakov = create_rw_buffer(sizeof(hmc_complex) * global_work_size);

	// scratch buffers for gauge observable will be created on demand
	clmem_plaq_buf_glob = 0;
	clmem_tplaq_buf_glob = 0;
	clmem_splaq_buf_glob = 0;
	clmem_polyakov_buf_glob = 0;

	return HMC_SUCCESS;
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
	if(get_parameters()->get_use_smearing()==true){
		stout_smear = createKernel("stout_smear") << basic_opencl_code << "gaugemomentum.cl" << "stout_smear.cl";
	}
	
}

hmc_error Opencl::init(cl_device_type wanted_device_type, inputparameters* params)
{
	hmc_error err = init_basic(wanted_device_type, params);
	return err;
}

hmc_error Opencl::init_basic(cl_device_type wanted_device_type, inputparameters* params)
{
	//variables, initializing, ...
	set_parameters(params);
	hmc_error err;
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

        if( strstr( info, "cl_amd_fp64" ) != NULL ) device_double_extension="AMD";
        if( strstr( info, "cl_khr_fp64" ) != NULL ) device_double_extension="KHR";


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
#ifdef _PROFILING_	
	queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &clerr);
#else
	queue = clCreateCommandQueue(context, device, 0, &clerr);
#endif
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "... failed, aborting.";
		exit(HMC_OCLERROR);
	}

	//Create buffer
	err = this->fill_buffers();
	if( err )
		exit( HMC_OCLERROR );

	//Create kernels
	this->fill_kernels();

	//finish
	set_init_true();
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

	if(clReleaseKernel(plaquette) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(plaquette_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(clReleaseKernel(polyakov_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(get_parameters()->get_use_smearing()==true){
		if(clReleaseKernel(stout_smear) != CL_SUCCESS) exit(HMC_OCLERROR);
	}

	return HMC_SUCCESS;
}

hmc_error Opencl::clear_buffers()
{
	logger.trace() << "Clearing memory objects";

	if(clReleaseMemObject(clmem_gaugefield) != CL_SUCCESS) exit(HMC_OCLERROR);

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

hmc_error Opencl::copy_buffer_on_device(cl_mem in, cl_mem out, size_t size, usetimer* timer)
{
	(*timer).reset();
	int clerr = CL_SUCCESS;

	clerr = clEnqueueCopyBuffer(queue, in, out, 0, 0, size , 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... copying buffer on device failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	(*timer).add();
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

#ifdef _PROFILING_
	cl_event event;
	clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, &event);
	int done = clWaitForEvents(1, &event);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clWaitForEvents failed with errorcode " << done << "aborting...";
		exit (HMC_OCLERROR);
	}
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
// 	char timerName[7] = ("timer_");
// 	strcat(timerName, kernelName);
// 	//Problem: How to call the memory object?
// 	(this->timerName).add(get_kernel_exec_time(event));
	
	//Second Method: Nasty workaround
	//noop is used in case the kernel is not recognized
	usetimer *noop = NULL;
	noop = Opencl::get_timer(kernelName);
	if(noop == NULL) 
		logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
	else
 		(*get_timer(kernelName)).add(get_kernel_exec_time(event));
#else
	clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, local_work_size_p, 0, 0, NULL);
#endif
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
	cl_int clerr;
#ifdef _PROFILING_
	cl_event event;
	clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, &event);
	int done = clWaitForEvents(1, &event);
	if(clerr != CL_SUCCESS) {
		logger.fatal() << "clWaitForEvents failed with errorcode " << done << "aborting...";
		exit (HMC_OCLERROR);
	}
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
// 	char timerName[7] = ("timer_");
// 	strcat(timerName, kernelName);
// 	//Problem: How to call the memory object?
// 	(this->timerName).add(get_kernel_exec_time(event));
	
	//Second Method: Nasty workaround
	//noop is used in case the kernel is not recognized
	usetimer *noop = NULL;
// 	noop = Opencl::get_timer(kernelName);
	noop = get_timer(kernelName);
	if(noop == NULL) 
		logger.error() << "get_timer(" << kernelName << ") did not return a timer!";
	else
 		(*get_timer(kernelName)).add(get_kernel_exec_time(event));
#else
	clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
#endif	
	
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

hmc_error Opencl::gaugeobservables(cl_mem gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	cl_int clerr = CL_SUCCESS;

	// decide on work-sizes
	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;
	//CP: This has no effect yet!!
	char * kernelname = "dummy";
	Opencl::get_work_sizes(&local_work_size, &global_work_size, &num_groups, device_type, kernelname);

	logger.debug() <<"init scratch buffers if not already done";
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;

	if( clmem_plaq_buf_glob == 0 ) clmem_plaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_tplaq_buf_glob == 0 ) clmem_tplaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_splaq_buf_glob == 0 ) clmem_splaq_buf_glob = create_rw_buffer(global_buf_size_float);
	if( clmem_polyakov_buf_glob == 0 ) clmem_polyakov_buf_glob = create_rw_buffer(global_buf_size_complex);
	
	//measure plaquette

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

	//measure polyakovloop
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

	return HMC_SUCCESS;
}

TmpClKernel Opencl::createKernel(const char * const kernel_name)
{
	stringstream collect_options;
	this->fill_collect_options(&collect_options);
	return TmpClKernel(kernel_name, collect_options.str(), context, &device, 1);
}

hmc_error Opencl::stout_smear_device(const size_t ls, const size_t gs){
	
	return HMC_SUCCESS;
}

hmc_error Opencl::get_work_sizes(size_t * ls, size_t * gs, cl_uint * num_groups, cl_device_type dev_type, char * name){
	/// @todo use kernelname
	size_t local_work_size;
	int numthreads = 128;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		local_work_size = numthreads; /// @todo have local work size depend on kernel properties (and device? autotune?)
	else
		local_work_size = 1; // nothing else makes sense on CPU

	size_t global_work_size;
	if( dev_type == CL_DEVICE_TYPE_GPU )
		global_work_size = 4 * numthreads * max_compute_units; /// @todo autotune
	else
		global_work_size = max_compute_units;

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;
	
	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;
	
	return HMC_SUCCESS;
}

#ifdef _PROFILING_
usetimer* Opencl::get_timer(char * in){
	if (strcmp(in, "polyakov_reduction") == 0){
    return &(this->timer_polyakov_reduction);
	}
	if (strcmp(in, "polyakov") == 0){
    return &(this->timer_polyakov);
	}
	if (strcmp(in, "plaquette_reduction") == 0){
    return &(this->timer_plaquette_reduction);
	}
	if (strcmp(in, "plaquette") == 0){
    return &(this->timer_plaquette);
	}
	if (strcmp(in, "stout_smear") == 0){
    return &(this->timer_stout_smear);
	}
	//if the kernelname has not matched, return NULL
	else{
		return NULL;
	}
}

int Opencl::get_read_write_size(char * in, inputparameters * parameters){
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	int R = (*parameters).get_mat_size();
	int S;
	if((*parameters).get_use_eo() == 1)
	  S = EOPREC_SPINORFIELDSIZE;
	else
	  S = SPINORFIELDSIZE;
	//this is the same as in the function above
	if (strcmp(in, "polyakov") == 0){
    return VOL4D*D*R + 1;
	}
	if (strcmp(in, "polyakov_reduction") == 0){
		//this is not right, since one does not know bufelements now
		//return (Bufel + 1) *2
    return NUMTHREADS;
	}
	if (strcmp(in, "plaquette") == 0){
    return 48*VOL4D *D*R + 1;
	}
	if (strcmp(in, "plaquette_reduction") == 0){
		//this is not right, since one does not know bufelements now
		//return (Bufel + 1) *2
    return NUMTHREADS;	
	}
	if (strcmp(in, "stout_smear") == 0){
    return 1000000000000000000000;
	}
	return 0;
}

void Opencl::print_profiling(std::string filename, char * kernelName, uint64_t time_total, int calls_total, int read_write_size){
	hmc_float bandwidth = 0.;
	uint64_t avg_time = 0.;
	uint64_t avg_time_site = 0.;
	//check if kernel has been called at all
	if(calls_total != 0 && time_total != 0){
		avg_time = (uint64_t) ( ( (float) time_total ) / ((float) calls_total) );
		avg_time_site =(uint64_t) ( ( (float) time_total ) / ((float) (calls_total*VOL4D)) );
		//Bandwidth in GB/s: 1e-3 = 1e6 (museconds) * 1e-9 (GByte)
		bandwidth = (hmc_float) read_write_size / (hmc_float) time_total * (hmc_float) calls_total *1e-3;
	}
	float mega = 1024*1024;
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) exit(HMC_FILEERROR);
	out.width(8);
	out.precision(15);
	//to look like that
	/*
	logger.trace() << "*******************************************************************";
	logger.trace() << "Fermion\t"<< setfill(' ') << setw(16)<< "BW[GB/s]\t" << setfill(' ') << setw(18) << "Re/Wr[MByte]\t" << setfill(' ') << setw(6)  << "Calls\t" << setfill(' ') << setw(10)  << "Time[mus]";
	*/
	out << kernelName << "\t" << time_total << "\t" << calls_total << "\t" << avg_time << "\t" << avg_time_site << "\t" << bandwidth << "\t" << (float) read_write_size/mega << std::endl;
	out.close();
	return;
}

void Opencl::print_profiling(std::string filename){
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

