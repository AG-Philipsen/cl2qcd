#include "opencl.h"

#include <algorithm>

using namespace std;

hmc_error opencl::init(cl_device_type wanted_device_type, const size_t local_work_size, const size_t global_work_size, usetimer* timer, inputparameters* parameters)
{

//give a list of all kernel-files
//!!CP: LZ should update this
	cl_kernels_file.push_back("opencl_header.cl");
	cl_kernels_file.push_back("opencl_geometry.cl");
	cl_kernels_file.push_back("opencl_random.cl");
	cl_kernels_file.push_back("opencl_operations_complex.cl");
	cl_kernels_file.push_back("opencl_operations_gaugefield.cl");
	cl_kernels_file.push_back("opencl_update_heatbath.cl");
	cl_kernels_file.push_back("opencl_gaugeobservables.cl");
#ifdef _FERMIONS_
	cl_kernels_file.push_back("opencl_operations_spinor.cl");
	cl_kernels_file.push_back("opencl_operations_spinorfield.cl");
	cl_kernels_file.push_back("opencl_operations_fermionmatrix.cl");
	cl_kernels_file.push_back("opencl_fermionobservables.cl");
#endif
#ifdef _TESTING_
	cl_kernels_file.push_back("opencl_testing.cl");
#endif

	cl_int clerr = CL_SUCCESS;

	timer->reset();
	cout << "OpenCL being initialized..." << endl;
	cl_uint num_platforms;
	cl_platform_id platform;
	//LZ: for now, stick to one platform without any further checks...
	clerr = clGetPlatformIDs(1, &platform, &num_platforms);
	if(clerr != CL_SUCCESS) {
		cout << "clGetPlatformIDs failed..." << endl;
		exit(HMC_OCLERROR);
	}

	char info[512];
	if(clGetPlatformInfo(platform, CL_PLATFORM_NAME, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\tCL_PLATFORM_NAME:     " << info << endl;
	if(clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\tCL_PLATFORM_VENDOR:   " << info << endl;
	if(clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\tCL_PLATFORM_VERSION:  " << info << endl;
	cout << endl;

	cl_uint num_devices;
	cl_device_id device;
	clerr = clGetDeviceIDs(platform, wanted_device_type, 0, NULL, &num_devices);
	if(num_devices == 1) {
		cout << "\t" << num_devices << " device of wanted type has been found." << endl;
	} else {
		cout << "\t" << num_devices << " devices of wanted type have been found. Choosing device number " << 0 << "." << endl;
	}
	clerr = clGetDeviceIDs(platform, wanted_device_type, 1, &device, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clGetDeviceIDs failed..." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "\tDevice information: " << endl;
	if(clGetDeviceInfo(device, CL_DEVICE_NAME, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\t\tCL_DEVICE_NAME:    " << info << endl;
	if(clGetDeviceInfo(device, CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\t\tCL_DEVICE_VENDOR:  " << info << endl;
	cl_device_type type;
	if(clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(cl_device_type), &type, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	if(type == CL_DEVICE_TYPE_CPU) cout << "\t\tCL_DEVICE_TYPE:    CPU" << endl;
	if(type == CL_DEVICE_TYPE_GPU) cout << "\t\tCL_DEVICE_TYPE:    GPU" << endl;
	if(type == CL_DEVICE_TYPE_ACCELERATOR) cout << "\t\tCL_DEVICE_TYPE:    ACCELERATOR" << endl;
	if(type != CL_DEVICE_TYPE_CPU && type != CL_DEVICE_TYPE_GPU && type != CL_DEVICE_TYPE_ACCELERATOR) {
		cout << "unexpected CL_DEVICE_TYPE..." << endl;
		exit(HMC_OCLERROR);
	}
	if(clGetDeviceInfo(device, CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\t\tCL_DEVICE_VERSION: " << info << endl;
	if(clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, 512 * sizeof(char), info, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
	cout << "\t\tCL_DEVICE_EXTENSIONS: " << info << endl;

	// figure out the number of "cores"
	if(clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL) != CL_SUCCESS) exit(HMC_OCLERROR);

	cout << "Create context..." << endl;
	context = clCreateContext(0, 1, &device, 0, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Create command queue..." << endl;
	queue = clCreateCommandQueue(context, device, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	// create array to point to contents of the different source files
	char ** sources = new char *[ cl_kernels_file.size() ];
	size_t * source_sizes = new size_t[ cl_kernels_file.size() ];

	string sourcecode;
	for(size_t n = 0; n < cl_kernels_file.size(); n++) {
		stringstream tmp;
		tmp << SOURCEDIR << '/' << cl_kernels_file[n];
		cout << "Read kernel source from file: " << tmp.str() << endl;

		fstream kernelsfile;
		kernelsfile.open(tmp.str().c_str());
		if(!kernelsfile.is_open()) {
			cerr << "Could not open file. Aborting..." << endl;
			exit(HMC_FILEERROR);
		}

		kernelsfile.seekg(0, ios::end);
		source_sizes[n] = kernelsfile.tellg();
		kernelsfile.seekg(0, ios::beg);

		sources[n] = new char[source_sizes[n]];

		kernelsfile.read( sources[n], source_sizes[n] );

		kernelsfile.close();
	}

	cout << "Create program..." << endl;
	clprogram = clCreateProgramWithSource(context, cl_kernels_file.size() , (const char**) sources, source_sizes, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Build program..." << endl;
	stringstream collect_options;
	collect_options << "-D_INKERNEL_ -DNSPACE=" << NSPACE << " -DNTIME=" << NTIME << " -DVOLSPACE=" << VOLSPACE << " -DSPINORSIZE=" << SPINORSIZE << " -DHALFSPINORSIZE=" << HALFSPINORSIZE << " -DSPINORFIELDSIZE=" << SPINORFIELDSIZE << " -DEOPREC_SPINORFIELDSIZE=" << EOPREC_SPINORFIELDSIZE;

	//CP: these have to match those in the cmake file
#ifdef _RECONSTRUCT_TWELVE_
	collect_options << " -D_RECONSTRUCT_TWELVE_";
#endif
#ifdef _USEDOUBLEPREC_
	collect_options << " -D_USEDOUBLEPREC_";
#endif
#ifdef _FERMIONS_
	collect_options << " -D_FERMIONS_";
#ifdef _TWISTEDMASS_
	collect_options << " -D_TWISTEDMASS_";
#endif
#ifdef _CLOVER_
	collect_options << " -D_CLOVER_";
#endif
#endif //_FERMIONS_
#ifdef _USEGPU_
	collect_options << " -D_USEGPU_";
#endif
#ifdef _PERFORM_BENCHMARKS_
	collect_options << " -D_PERFORM_BENCHMARKS_";
#endif
#ifdef _CP_REAL_
	collect_options << " -D_CP_REAL_";
#endif
#ifdef _CP_IMAG_
	collect_options << " -D_CP_IMAG_";
#endif
#ifdef _NPBC_T_
	collect_options << " -D_NPBC_T_";
#endif
#ifdef _NPBC_S_
	collect_options << " -D_NPBC_S_";
#endif

#ifdef _FERMIONS_
	//CP: give kappa and its negative value
	hmc_float kappa_tmp = parameters->get_kappa();
	collect_options << " -DKAPPA=" << kappa_tmp;
	collect_options << " -DMKAPPA=" << -kappa_tmp;
#ifdef _TWISTEDMASS_
	hmc_float mu_tmp = parameters->get_mu();
	collect_options << " -DMU=" << mu_tmp;
#endif
#ifdef _CLOVER_
	hmc_float csw_tmp = parameters->get_csw();
	collect_options << " -DCSW=" << csw_tmp;
#endif
#endif //_FERMIONS_

	collect_options << " -I" << SOURCEDIR;
	string buildoptions = collect_options.str();
	cout << "\tbuild options:";
	cout << "\t" << buildoptions << endl;
	clerr = clBuildProgram(clprogram, 1, &device, buildoptions.c_str(), 0, 0);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, but look at BuildLog and abort then." << endl;
	}

	cout << "finished building program" << endl;
	size_t logSize;
	clerr |= clGetProgramBuildInfo(clprogram, device, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
	if(logSize > 1) { // 0-terminated -> always at least one byte
		cout << "Build Log:" << endl;
		char* log = new char[logSize];
		clerr |= clGetProgramBuildInfo(clprogram, device, CL_PROGRAM_BUILD_LOG, logSize, log, NULL);
		cout << log << endl << endl;
		delete [] log;
	}
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;

		// dump program source
		size_t sourceSize;
		clerr = clGetProgramInfo(clprogram, CL_PROGRAM_SOURCE, 0, NULL, &sourceSize);
		if(!clerr && sourceSize > 1) { // 0-terminated -> always at least one byte
			char* source = new char[sourceSize];
			clerr = clGetProgramInfo(clprogram, CL_PROGRAM_SOURCE, sourceSize, source, &sourceSize);
			if(!clerr) {
				char const * const FILENAME = "broken_source.cl";
				ofstream srcFile(FILENAME);
				srcFile << source;
				srcFile.close();
				cout << "Dumped broken source to " << FILENAME << endl;
			}
			delete[] source;
		}

		exit(HMC_OCLERROR);
	}

	cout << "Create buffer for gaugefield..." << endl;
	clmem_gaugefield = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_gaugefield), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Create buffer for random numbers..." << endl;
	clmem_rndarray = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_rndarray), 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Create buffer for gaugeobservables..." << endl;
	clmem_plaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clmem_splaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clmem_tplaq = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clmem_polyakov = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex) * global_work_size, 0, &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	// scratch buffers for gauge observable will be created on demand
	clmem_plaq_buf_glob = 0;
	clmem_tplaq_buf_glob = 0;
	clmem_splaq_buf_glob = 0;
	clmem_polyakov_buf_glob = 0;

	cout << "Create heatbath kernels..." << endl;
	heatbath_even = clCreateKernel(clprogram, "heatbath_even", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	heatbath_odd = clCreateKernel(clprogram, "heatbath_odd", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	overrelax_even = clCreateKernel(clprogram, "overrelax_even", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	overrelax_odd = clCreateKernel(clprogram, "overrelax_odd", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	cout << "Create gaugeobservables kernels..." << endl;
	plaquette = clCreateKernel(clprogram, "plaquette", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	polyakov = clCreateKernel(clprogram, "polyakov", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	plaquette_reduction = clCreateKernel(clprogram, "plaquette_reduction", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	polyakov_reduction = clCreateKernel(clprogram, "polyakov_reduction", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	isinit = 1;

	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_gaugefield_to_device(hmc_gaugefield* gaugefield, usetimer* timer)
{
//   cout<<"Copy gaugefield to device..."<<endl;
	timer->reset();
	hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

	copy_to_ocl_format(host_gaugefield, gaugefield);

	int clerr = clEnqueueWriteBuffer(queue, clmem_gaugefield, CL_TRUE, 0, sizeof(hmc_gaugefield), host_gaugefield, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	free(host_gaugefield);

	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_rndarray_to_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Copy randomarray to device..."<<endl;
	timer->reset();

	int clerr = clEnqueueWriteBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::get_gaugefield_from_device(hmc_gaugefield* gaugefield, usetimer* timer)
{
//   cout<<"Get gaugefield from device..."<<endl;
	timer->reset();
	hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

	int clerr = clEnqueueReadBuffer(queue, clmem_gaugefield, CL_TRUE, 0, sizeof(hmc_gaugefield), host_gaugefield, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		cout << "errorcode :" << clerr << endl;
		exit(HMC_OCLERROR);
	}

	copy_from_ocl_format(gaugefield, host_gaugefield);

	free(host_gaugefield);

	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer)
{
//   cout<<"Get randomarray from device..."<<endl;
	timer->reset();

	int clerr = clEnqueueReadBuffer(queue, clmem_rndarray, CL_TRUE, 0, sizeof(hmc_rndarray), rndarray, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}

	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::run_heatbath(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;
	timer->reset();

#ifdef _USE_GPU_
	size_t global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
#else
	size_t global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);
#endif

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 at heatbath_even failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 at heatbath_even failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 at heatbath_even failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			cout << "clSetKernelArg2 at heatbath_even failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_even, global_work_size);
	}

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 at heatbath_odd failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 at heatbath_odd failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(heatbath_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 at heatbath_odd failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			cout << "clSetKernelArg2 at heatbath_odd failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
		enqueueKernel(heatbath_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;

}

hmc_error opencl::run_overrelax(const hmc_float beta, usetimer * const timer)
{
	cl_int clerr = CL_SUCCESS;

	timer->reset();

#ifdef _USE_GPU_
	size_t global_work_size = min(VOLSPACE * NTIME / 2, NUMRNDSTATES);
#else
	size_t global_work_size = min(max_compute_units, (cl_uint) NUMRNDSTATES);
#endif

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_even, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			cout << "clSetKernelArg4 failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_even, global_work_size);
	}

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), &clmem_gaugefield);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg6 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(overrelax_odd, 3, sizeof(cl_mem), &clmem_rndarray);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg7 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	for(int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 2, sizeof(int), &i);
		if(clerr != CL_SUCCESS) {
			cout << "clSetKernelArg8 failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
		enqueueKernel(overrelax_odd, global_work_size);
	}
	clFinish(queue);
	timer->add();
	return HMC_SUCCESS;
}

hmc_error opencl::gaugeobservables(hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, usetimer* timer1, usetimer * timer2)
{
	cl_int clerr = CL_SUCCESS;

	// decide on work-sizes
#ifdef _USE_GPU_
	const size_t local_work_size = NUM_THREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t local_work_size = 1; // nothing else makes sens on CPU
#endif

#ifdef _USE_GPU_
	size_t global_work_size = 4 * NUM_THREADS * max_compute_units; /// @todo autotune
#else
	size_t global_work_size = max_compute_units;
#endif

	const cl_uint num_groups = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups;

	// init scratch buffers if not already done
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;

	if( clmem_plaq_buf_glob == 0 ) {
		clmem_plaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating clmem_plaq_buf_glob failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_tplaq_buf_glob == 0 ) {
		clmem_tplaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating tclmem_plaq_buf_glob failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_splaq_buf_glob == 0 ) {
		clmem_splaq_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_float, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating sclmem_plaq_buf_glob failed, aborting..." << endl;
			exit(HMC_OCLERROR);
		}
	}
	if( clmem_polyakov_buf_glob == 0 ) {
		clmem_polyakov_buf_glob = clCreateBuffer(context, CL_MEM_READ_WRITE, global_buf_size_complex, 0, &clerr);
		if(clerr != CL_SUCCESS) {
			cout << "creating clmem_polyakov_buf_glob failed, aborting..." << endl;
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
		cout << "clSetKernelArg0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg6 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette, global_work_size, local_work_size);

	// run second part of plaquette reduction

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), &clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), &clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), &clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), &clmem_plaq);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg3 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), &clmem_tplaq);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg4 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), &clmem_splaq);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg5 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg6 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(plaquette_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_plaq, CL_FALSE, 0, sizeof(hmc_float), &plaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_tplaq, CL_FALSE, 0, sizeof(hmc_float), &tplaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clEnqueueReadBuffer(queue, clmem_splaq, CL_FALSE, 0, sizeof(hmc_float), &splaq, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
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
		cout << "clSetKernelArg0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov, global_work_size, local_work_size);

	// second part of polyakov reduction

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), &clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg0 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), &clmem_polyakov);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg1 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) {
		cout << "clSetKernelArg2 failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
	enqueueKernel(polyakov_reduction, 1, 1);

	//read out values
	clerr = clEnqueueReadBuffer(queue, clmem_polyakov, CL_FALSE, 0, sizeof(hmc_complex), &pol, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "... failed, aborting." << endl;
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


hmc_error opencl::finalize()
{
	if(isinit == 1) {
		if(clFlush(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clFinish(queue) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(heatbath_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(heatbath_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(overrelax_even) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(overrelax_odd) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseKernel(plaquette) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(polyakov) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(plaquette_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseKernel(polyakov_reduction) != CL_SUCCESS) exit(HMC_OCLERROR);

		if(clReleaseProgram(clprogram) != CL_SUCCESS) exit(HMC_OCLERROR);

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

		if(clReleaseCommandQueue(queue) != CL_SUCCESS) exit(HMC_OCLERROR);
		if(clReleaseContext(context) != CL_SUCCESS) exit(HMC_OCLERROR);

		isinit = 0;
	}
	return HMC_SUCCESS;
}

void opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size)
{
#ifdef _USE_GPU_
	const size_t local_work_size = NUM_THREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t local_work_size = 1; // nothing else makes sens on CPU
#endif

	cl_int clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clEnqueueNDRangeKernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
}

void opencl::enqueueKernel(const cl_kernel kernel, const size_t global_work_size, const size_t local_work_size)
{
	cl_int clerr = clEnqueueNDRangeKernel(queue, kernel, 1, 0, &global_work_size, &local_work_size, 0, 0, NULL);
	if(clerr != CL_SUCCESS) {
		cout << "clEnqueueNDRangeKernel failed, aborting..." << endl;
		exit(HMC_OCLERROR);
	}
}

//CP: extra file for the long testing functions
#include "opencl_testing.h"
//CP: extra file for the long fermion functions
#ifdef _FERMIONS_
#include "opencl_fermions.h"
#endif //_FERMIONS_
