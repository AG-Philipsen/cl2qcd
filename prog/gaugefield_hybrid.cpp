#include "gaugefield_hybrid.h"

void Gaugefield_hybrid::init(int numtasks, cl_device_type primary_device_type, inputparameters* input_parameters)
{

	logger.trace() << "Initialize gaugefield";

	//how many tasks (devices with different purpose) do we need:
	set_num_tasks(numtasks);
	//LZ: for now assume that there is only one device per task

	//input parameters
	set_parameters(input_parameters);

	//allocate memory for private gaugefield on host and initialize (cold start, read in, future: hot start)
	sgf = new Matrixsu3[get_num_gaugefield_elems()];
	init_gaugefield();

	init_devicetypearray(primary_device_type);
	init_opencl();

	this->init_tasks();

	//this has to be done anyways...
	copy_gaugefield_to_all_tasks();
}

void Gaugefield_hybrid::init_devicetypearray(cl_device_type primary_device_type)
{
	//init devicetype array

	//LZ: Note that num_task_types and the input parameter num_dev seem not to fit to each other. However, consider following scenario:
	//    num_task_types is given to this class and can potentially control whether certain additional tasks are performed or not
	//    num_dev controls how many devices should be used/are available (so usually CPU and GPU, but possibly also CPU, GPU from different nodes)
	//    here, the different devices can be assigned automatically to the tasks
	//CP: At the moment, get_num_dev is simply not used, altough it may be in the future...
	//LZ: first application of num_dev: if num_dev==1 then we want to use devices of one type (primary_device_type) only

	devicetypes = new cl_device_type[get_num_tasks()];
	if(get_num_tasks() == 1) {
		devicetypes[0] = primary_device_type;
	} else if (get_num_tasks() == 2) {
		devicetypes[0] = primary_device_type;
		if(primary_device_type == CL_DEVICE_TYPE_GPU) {
			if(get_parameters()->get_num_dev() == 1)
				devicetypes[1] = CL_DEVICE_TYPE_GPU;
			else
				devicetypes[1] = CL_DEVICE_TYPE_CPU;
		} else {
			if(get_parameters()->get_num_dev() == 1)
				devicetypes[1] = CL_DEVICE_TYPE_CPU;
			else
				devicetypes[1] = CL_DEVICE_TYPE_GPU;
		}
	} else {
		throw Print_Error_Message("3 or more tasks not yet implemented.");
	}

	logger.debug() << "Wish list for device types:" ;
	for(int n = 0; n < get_num_tasks(); n++) {
		switch(devicetypes[n]) {
			case CL_DEVICE_TYPE_CPU :
				logger.debug() << "CL_DEVICE_TYPE_CPU" ;
				break;
			case CL_DEVICE_TYPE_GPU :
				logger.debug() << "CL_DEVICE_TYPE_GPU" ;
				break;
		}
	}
}

void Gaugefield_hybrid::init_opencl()
{
	cl_int clerr = CL_SUCCESS;

	// in debug scenarios make the compiler dump the compile results
	if( logger.beDebug() ) {
		setenv("GPU_DUMP_DEVICE_KERNEL", "3", 0); // can be overriden from outside
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


	cl_uint num_devices_gpu;
	cl_uint num_devices_cpu;
	clerr = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices_gpu);
	clerr = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &num_devices_cpu);
	logger.info() << "\tFound " << num_devices_gpu << " GPU(s) and " << num_devices_cpu << " CPU(s).";

	//LZ: begin debug
	//    num_devices_gpu = 0;
	//  num_devices_cpu = 1;
	//  logger.info() << "Found " << num_devices_gpu << " GPU(s) and " << num_devices_cpu << " CPU(s).";
	// end debug

	set_num_devices(num_devices_gpu + num_devices_cpu);

	//Map tasks and devices
	if(get_num_devices() < 1) throw Print_Error_Message("No devices found.");
	if( num_devices_gpu == 0 ) logger.warn() << "No GPU found.";
	if( num_devices_cpu == 0 ) logger.warn() << "No CPU found.";

	switch (get_num_tasks() ) {

		case 1 :
			if(devicetypes[0] == CL_DEVICE_TYPE_GPU && num_devices_gpu < 1) {
				logger.warn() << "Program demanded GPU device but there is only a CPU device available. Try that." ;
				devicetypes[0] = CL_DEVICE_TYPE_CPU;
			}
			if(devicetypes[0] == CL_DEVICE_TYPE_CPU && num_devices_cpu < 1) {
				logger.warn() << "Program demanded CPU device but there is only a GPU device available. Try that." ;
				devicetypes[0] = CL_DEVICE_TYPE_GPU;
			}
			break;

		case 2 :
			if( get_num_devices() == 1 ) {
				logger.warn() << "You wanted to have two devices, but only one has been found!";
				if( num_devices_gpu == 1 ) {
					devicetypes[0] = CL_DEVICE_TYPE_GPU;
					devicetypes[1] = CL_DEVICE_TYPE_GPU;
				}
				if( num_devices_cpu == 1 ) {
					devicetypes[0] = CL_DEVICE_TYPE_CPU;
					devicetypes[1] = CL_DEVICE_TYPE_CPU;
				}
			}
			break;

		default:
			logger.warn() << "Set of devices could not be mapped properly. Try what happens..." ;
	}

	int len = min( get_num_tasks(), get_num_devices() );
	queue                   = new cl_command_queue [get_num_tasks()];
	devices                 = new cl_device_id     [len];
	device_double_extension = new string  [len];
	max_compute_units       = new cl_uint [len];

	for(int ntask = 0; ntask < len; ntask++) {
		logger.info() << "\tInitialize device #" << ntask << ":";
		init_devices(ntask);
	}

	//Initilize context
	logger.trace() << "Create context...";
	context = clCreateContext(0, len, devices, 0, 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateContext", __FILE__, __LINE__);

	//now we need a mapping between devices and tasks for the case of fewer devices than tasks
	device_id_for_task = new int [len];
	for(int ntask = 0 ; ntask < get_num_tasks() ; ntask++) {
		device_id_for_task[ntask] = ntask;
	}
	if(get_num_devices() < get_num_tasks()) {
		switch (get_num_tasks()) {
			case 2:
				device_id_for_task[1] = 0;
				break;
			default:
				throw Print_Error_Message("Fewer devices than tasks but no proper mapping available.");
		}
	}

	//Initilize queues, one per task
	// Note that it might be advantageous to combine tasks on the same device into the same queue, i.e. to have only one queue per device even for more devices
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		logger.trace() << "Create command queue for task #" << ntask << "...";

#ifdef _PROFILING_
		queue[ntask] = clCreateCommandQueue(context, get_device_for_task(ntask), CL_QUEUE_PROFILING_ENABLE, &clerr);
#else
		queue[ntask] = clCreateCommandQueue(context, get_device_for_task(ntask), 0, &clerr);
#endif
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateCommandQueue", __FILE__, __LINE__);
	}


	//context-wide buffers
	logger.trace() << "Creating gaugefield buffer...";
	clmem_gaugefield = clCreateBuffer(context, CL_MEM_READ_WRITE, get_num_gaugefield_elems() * sizeof(Matrixsu3), 0, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);
}

void Gaugefield_hybrid::init_devices(int ndev)
{
	cl_int clerr = CL_SUCCESS;

	char info[512];

	clerr = clGetDeviceIDs(platform, get_device_type(ndev), 1, &devices[ndev], NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceIDs", __FILE__, __LINE__);

	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_NAME, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.info() << "\t\t\tCL_DEVICE_NAME:    " << info;
	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.debug() << "\t\t\tCL_DEVICE_VENDOR:  " << info;
	cl_device_type devtype;
	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_TYPE, sizeof(cl_device_type), &devtype, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	if(devtype == CL_DEVICE_TYPE_CPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    CPU";
	if(devtype == CL_DEVICE_TYPE_GPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    GPU";
	if(devtype == CL_DEVICE_TYPE_ACCELERATOR) logger.info() << "\t\t\tCL_DEVICE_TYPE:    ACCELERATOR";
	if(devtype != CL_DEVICE_TYPE_CPU && devtype != CL_DEVICE_TYPE_GPU && devtype != CL_DEVICE_TYPE_ACCELERATOR)
		throw Print_Error_Message("Unexpected CL_DEVICE_TYPE...", __FILE__, __LINE__);
	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.debug() << "\t\t\tCL_DEVICE_VERSION: " << info;
	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_EXTENSIONS, 512 * sizeof(char), info, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);
	logger.debug() << "\t\t\tCL_DEVICE_EXTENSIONS: " << info;

	if( strstr( info, "cl_amd_fp64" ) != NULL ) device_double_extension[ndev] = "AMD";
	if( strstr( info, "cl_khr_fp64" ) != NULL ) device_double_extension[ndev] = "KHR";

	// figure out the number of "cores"
	clerr = clGetDeviceInfo(devices[ndev], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units[ndev], NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetDeviceInfo", __FILE__, __LINE__);;
	logger.debug() << "\t\t\tCL_DEVICE_MAX_COMPUTE_UNITS: " << max_compute_units[ndev];
}



void Gaugefield_hybrid::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
		opencl_modules[ntask] = new Opencl_Module[1];
		opencl_modules[ntask]->init(queue[ntask], &clmem_gaugefield, get_parameters(), max_compute_units[ntask], get_double_ext(ntask));
	}
}


void Gaugefield_hybrid::finalize()
{
	this->finalize_opencl();
	this->delete_variables();
}

void Gaugefield_hybrid::delete_variables()
{
	delete [] device_id_for_task;

	delete [] sgf;

	delete [] devices;
	delete [] queue;
	delete [] device_double_extension;
	delete [] max_compute_units;

	delete [] devicetypes;

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		delete [] opencl_modules[ntask];
	}
	delete [] opencl_modules;
}

void Gaugefield_hybrid::finalize_opencl()
{
	cl_int clerr = CL_SUCCESS;

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		clerr = clFlush(queue[ntask]);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFlush", __FILE__, __LINE__);
		clerr = clFinish(queue[ntask]);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);
	}

	//  tasks->clear_kernels();
	//  tasks->clear_buffers();

	clerr = clReleaseMemObject(clmem_gaugefield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		clerr = clReleaseCommandQueue(queue[ntask]);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseCommandQueue", __FILE__, __LINE__);
	}

	clerr = clReleaseContext(context);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseContext", __FILE__, __LINE__);
}

void Gaugefield_hybrid::init_gaugefield()
{
	if((get_parameters())->get_startcondition() == START_FROM_SOURCE) {
		sourcefileparameters parameters_source;

		//hmc_gaugefield for filetransfer, initialize here, because otherwise it is not needed
 		hmc_complex* gftmp = new hmc_complex[get_num_hmc_gaugefield_elems()];
		//tmp gauge field
		hmc_float * gaugefield_tmp;
		gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float) * NDIM * NC * NC * parameters->get_nt() * parameters->get_volspace());
		parameters_source.readsourcefile(&(get_parameters()->sourcefile)[0], get_parameters()->get_prec(), &gaugefield_tmp);
		copy_gaugefield_from_ildg_format(gftmp, gaugefield_tmp, parameters_source.num_entries_source, parameters);
 		copy_gaugefield_to_s_gaugefield (get_sgf(), gftmp);
		free(gaugefield_tmp);
 		delete[] gftmp;
	}
	if(get_parameters()->get_startcondition() == COLD_START) {
		set_gaugefield_cold(get_sgf());
	}
	if(get_parameters()->get_startcondition() == HOT_START) {
		set_gaugefield_hot(get_sgf());
	}
}

void Gaugefield_hybrid::set_gaugefield_cold(Matrixsu3 * field)
{
	for(int t = 0; t < parameters->get_nt(); t++) {
		for(int n = 0; n < parameters->get_volspace(); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const Matrixsu3 tmp = unit_matrixsu3();
				set_to_gaugefield(field, mu, n, t, tmp);
			}
		}
	}
}

//Implement this
void Gaugefield_hybrid::set_gaugefield_hot(Matrixsu3 *)
{
	throw Print_Error_Message("Hot start not yet implemented.", __FILE__, __LINE__);
}

void Gaugefield_hybrid::copy_gaugefield_to_all_tasks()
{
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		copy_gaugefield_to_task(ntask);
	}
}

void Gaugefield_hybrid::copy_gaugefield_to_task(int ntask)
{
	if(ntask < 0 || ntask > get_num_tasks() ) {
		logger.warn() << "Index out of range, copy_gaugefield_to_device does nothing.";
		return;
	}

	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(get_num_gaugefield_elems() * sizeof(ocl_s_gaugefield));

	copy_to_ocl_format(host_gaugefield, get_sgf(), parameters);

	cl_int clerr = clEnqueueWriteBuffer(queue[ntask], clmem_gaugefield, CL_TRUE, 0, get_num_gaugefield_elems() * sizeof(ocl_s_gaugefield), host_gaugefield, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	free(host_gaugefield);
}

void Gaugefield_hybrid::synchronize(int ntask_reference)
{
	if(ntask_reference < 0 || ntask_reference > get_num_tasks() ) {
		logger.warn() << "Index out of range, synchronize_gaugefield does nothing.";
		return;
	}
	clFinish(queue[ntask_reference]);
	copy_gaugefield_from_task(ntask_reference);

	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		clFinish(queue[ntask]);
		if(ntask != ntask_reference) copy_gaugefield_to_task(ntask);
	}
}

void Gaugefield_hybrid::copy_gaugefield_from_task(int ntask)
{
	if(ntask < 0 || ntask > get_num_tasks() ) {
		logger.warn() << "Index out of range, copy_gaugefield_from_device does nothing.";
		return;
	}

	ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(get_num_gaugefield_elems() * sizeof(ocl_s_gaugefield));

	cl_int clerr = clEnqueueReadBuffer(queue[ntask], clmem_gaugefield, CL_TRUE, 0, get_num_gaugefield_elems() * sizeof(ocl_s_gaugefield), host_gaugefield, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	copy_from_ocl_format(get_sgf(), host_gaugefield, parameters);

	free(host_gaugefield);
}

cl_device_id Gaugefield_hybrid::get_device_for_task(int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("task index out of range", __FILE__, __LINE__);
	return devices[device_id_for_task[ntask]];
}

cl_mem* Gaugefield_hybrid::get_clmem_gaugefield()
{
	return &clmem_gaugefield;
}

void Gaugefield_hybrid::set_num_tasks (int num)
{
	num_tasks = num;
}

int Gaugefield_hybrid::get_num_tasks ()
{
	return num_tasks;
}

int Gaugefield_hybrid::get_max_compute_units(int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("rndarray index out of range", __FILE__, __LINE__);
	return max_compute_units[device_id_for_task[ntask]];
}

string Gaugefield_hybrid::get_double_ext(int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("rndarray index out of range", __FILE__, __LINE__);
	return device_double_extension[device_id_for_task[ntask]];
}

inputparameters * Gaugefield_hybrid::get_parameters ()
{
	return  parameters;
}

void Gaugefield_hybrid::set_parameters (inputparameters * parameters_val)
{
	parameters = parameters_val;
}


Matrixsu3 * Gaugefield_hybrid::get_sgf ()
{
	return sgf;
}

void Gaugefield_hybrid::set_sgf (Matrixsu3 * sgf_val)
{
	sgf = sgf_val;
}

void Gaugefield_hybrid::set_num_devices(int num)
{
	num_devices = num;
}

int Gaugefield_hybrid::get_num_devices()
{
	return num_devices;
}

cl_device_type Gaugefield_hybrid::get_device_type(int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	return devicetypes[ntask];
}

void Gaugefield_hybrid::copy_gaugefield_to_s_gaugefield (Matrixsu3 * sgfo, hmc_complex * gf)
{
	const size_t NTIME = parameters->get_nt();
	for (int d = 0; d < NDIM; d++) {
		for (int n = 0; n < parameters->get_volspace(); n++) {
			for (size_t t = 0; t < NTIME; t++) {
				hmc_su3matrix srcElem;
				get_su3matrix(&srcElem, gf, n, t, d, parameters);
				Matrixsu3 destElem;
				destElem.e00 = srcElem[0][0];
				destElem.e01 = srcElem[0][1];
				destElem.e02 = srcElem[0][2];
				destElem.e10 = srcElem[1][0];
				destElem.e11 = srcElem[1][1];
				destElem.e12 = srcElem[1][2];
				destElem.e20 = srcElem[2][0];
				destElem.e21 = srcElem[2][1];
				destElem.e22 = srcElem[2][2];
				set_to_gaugefield(sgfo, d, n, t, destElem);
			}
		}
	}
}

void Gaugefield_hybrid::copy_s_gaugefield_to_gaugefield(hmc_complex * gf, Matrixsu3 * sgfo)
{
	const size_t NTIME = parameters->get_nt();
	for (int d = 0; d < NDIM; d++) {
		for (int n = 0; n < parameters->get_volspace(); n++) {
			for (size_t t = 0; t < NTIME; t++) {
				hmc_su3matrix destElem;
				Matrixsu3 srcElem = get_from_gaugefield(sgfo, d, n, t);
				destElem[0][0] = srcElem.e00;
				destElem[0][1] = srcElem.e01;
				destElem[0][2] = srcElem.e02;
				destElem[1][0] = srcElem.e10;
				destElem[1][1] = srcElem.e11;
				destElem[1][2] = srcElem.e12;
				destElem[2][0] = srcElem.e20;
				destElem[2][1] = srcElem.e21;
				destElem[2][2] = srcElem.e22;
				put_su3matrix(gf, &destElem, n, t, d, parameters);
			}
		}
	}
}


void Gaugefield_hybrid::save(int number)
{
	//LZ: generalize the following to larger numbers, if necessary...
	stringstream strnumber;
	strnumber.fill('0');
	strnumber.width(5);
	strnumber << right << number;
	stringstream outfilename;
	outfilename << "conf." << strnumber.str();
	string outputfile = outfilename.str();
	save(outputfile);
}


void Gaugefield_hybrid::save(string outputfile)
{
	const size_t NTIME = parameters->get_nt();
	const size_t gaugefield_buf_size = 2 * NC * NC * NDIM * parameters->get_volspace() * NTIME;
	hmc_float * gaugefield_buf = new hmc_float[gaugefield_buf_size];

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;

	hmc_complex* gftmp = new hmc_complex[get_num_hmc_gaugefield_elems()];
	copy_s_gaugefield_to_gaugefield(gftmp, get_sgf());
	copy_gaugefield_to_ildg_format(gaugefield_buf, gftmp, parameters);

	hmc_float plaq = plaquette();

	int number = 0;

	const size_t NSPACE = parameters->get_ns();

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NTIME, get_parameters()->get_prec(), number, plaq, get_parameters()->get_beta(), get_parameters()->get_kappa(), get_parameters()->get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	delete[] gaugefield_buf;
	delete[] gftmp;
}


hmc_float Gaugefield_hybrid::plaquette()
{
	hmc_float tdummy;
	hmc_float sdummy;
	//LZ: calculation of tdummy and sdummy is unnecessary for this function, chance to speed up a little bit...
	return plaquette(&tdummy, &sdummy);
}

void Gaugefield_hybrid::print_gaugeobservables(int iter)
{
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	hmc_complex pol = polyakov();
	logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
}

void Gaugefield_hybrid::print_gaugeobservables(int iter, std::string filename)
{
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	hmc_complex pol = polyakov();
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) throw File_Exception(filename);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	hmc_float plaq  = 0;
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_complex pol;
	cl_mem gf = *get_clmem_gaugefield();
	opencl_modules[ntask]->gaugeobservables(gf, &plaq, &tplaq, &splaq, &pol);
	logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask, std::string filename)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	hmc_float plaq  = 0;
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_complex pol;
	cl_mem gf = *get_clmem_gaugefield();
	opencl_modules[ntask]->gaugeobservables(gf, &plaq, &tplaq, &splaq, &pol);
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) throw File_Exception(filename);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
}

#ifdef _PROFILING_
void Gaugefield_hybrid::print_profiling(std::string filename)
{
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
		opencl_modules[ntask]->print_profiling(filename, ntask);
	}
}
#endif

size_t Gaugefield_hybrid::get_num_hmc_gaugefield_elems()
{
	return NC * NC * NDIM * parameters->get_volspace() * parameters->get_nt();
}

void Gaugefield_hybrid::set_to_gaugefield(Matrixsu3 * field, const size_t mu, const size_t x, const size_t t, const Matrixsu3 val)
{
	field[get_global_link_pos(mu, x, t, parameters)] = val;
}

Matrixsu3 Gaugefield_hybrid::get_from_gaugefield(const Matrixsu3 * field, const size_t mu, const size_t x, const size_t t) const
{
	return field[get_global_link_pos(mu, x, t, parameters)];
}

size_t Gaugefield_hybrid::get_num_gaugefield_elems() const
{
	return NDIM * parameters->get_volspace() * parameters->get_nt();
}

void Gaugefield_hybrid::copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, Matrixsu3* gaugefield, const inputparameters * const parameters)
{
	const size_t NSPACE = parameters->get_ns();
	const size_t NTIME = parameters->get_nt();
	for(size_t spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(size_t t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, parameters);
				host_gaugefield[index] = gaugefield[index];
			}
		}
	}
	return;
}

void Gaugefield_hybrid::copy_from_ocl_format(Matrixsu3* gaugefield, ocl_s_gaugefield* host_gaugefield, const inputparameters * const parameters)
{
	const size_t NSPACE = parameters->get_ns();
	const size_t NTIME = parameters->get_nt();
	for(size_t spacepos = 0; spacepos < NSPACE * NSPACE * NSPACE; spacepos++) {
		for(size_t t = 0; t < NTIME; t++) {
			for(int mu = 0; mu < NDIM; mu++) {
				const size_t index = get_global_link_pos(mu, spacepos, t, parameters);
				gaugefield[index] = host_gaugefield[index];
			}
		}
	}
	return;
}

void Gaugefield_hybrid::copy_gaugefield_from_ildg_format(hmc_complex * gaugefield, hmc_float * gaugefield_tmp, int check, const inputparameters * const parameters)
{
	//little check if arrays are big enough
	if (parameters->get_vol4d() *NDIM*NC*NC * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t VOLSPACE = parameters->get_volspace();
	const size_t NSPACE = parameters->get_ns();

	int cter = 0;
	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < parameters->get_nt(); t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (size_t i = 0; i < NSPACE; i++) {
			for (size_t j = 0; j < NSPACE; j++) {
				for (size_t k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in hmc_su3matrix tmp (NOTE: this is a complex array, not a struct!!)
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * VOLSPACE * NDIM;
						hmc_su3matrix tmp;
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers
								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								tmp[m][n].re = gaugefield_tmp[pos];
								tmp[m][n].im = gaugefield_tmp[pos + 1];
								cter++;
							}
						}
						put_su3matrix(gaugefield, &tmp, spacepos, t, (l + 1) % NDIM, parameters);
						/*
						//CP: This did not work because there is a non-trivial mapping going on in put_su3matrix
						//which needs to be cleared!!!
						//convert tmp to Matrixsu3 type and store it in the gaugefield
						Matrixsu3 tmp2;
						tmp2 = convert_hmc_matrixsu3_to_Matrixsu3(tmp);
						put_matrixsu3(gaugefield, tmp2, spacepos, t, (l + 1) % NDIM, parameters);
						*/
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	return;
}

void Gaugefield_hybrid::copy_gaugefield_to_ildg_format(hmc_float * dest, hmc_complex * source, const inputparameters * const parameters)
{

	int cter = 0;
	const size_t NSPACE = parameters->get_ns();

	//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
	for (int t = 0; t < parameters->get_nt(); t++) {
		//if the alg is known to be correct, the next three for-loops could also be changed to one
		for (size_t i = 0; i < NSPACE; i++) {
			for (size_t j = 0; j < NSPACE; j++) {
				for (size_t k = 0; k < NSPACE; k++) {
					for (int l = 0; l < NDIM; l++) {
						int spacepos = k + j * NSPACE + i * NSPACE * NSPACE;
						int globalpos = l + spacepos * NDIM + t * parameters->get_volspace() * NDIM;
						hmc_su3matrix tmp;
						get_su3matrix(&tmp, source, spacepos, t, (l + 1) % NDIM, parameters);
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
								//which is stored in one single array here
								//skip NC*NC*2 cmplx numbers
								int pos = 2 * n + 2 * m * NC + globalpos * NC * NC * 2;
								dest[pos]     = tmp[m][n].re;
								dest[pos + 1] = tmp[m][n].im;
								cter++;
							}
						}
					}
				}
			}
		}
	}

	return;
}

hmc_float Gaugefield_hybrid::plaquette(hmc_float* tplaq, hmc_float* splaq)
{
	hmc_float plaq = 0;
	*tplaq = 0;
	*splaq = 0;

	for(int t = 0; t < parameters->get_nt(); t++) {
		for(int n = 0; n < parameters->get_volspace(); n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				for(int nu = 0; nu < mu; nu++) {
					Matrixsu3 prod;
					prod = local_plaquette(get_sgf(), n, t, mu, nu, parameters );
					hmc_float tmpfloat = trace_Matrixsu3(prod).re;
					plaq += tmpfloat;
					if(mu == 0 || nu == 0) {
						*tplaq += tmpfloat;
					} else {
						*splaq += tmpfloat;
					}
				}
			}
		}
	}

	*tplaq /= static_cast<hmc_float>(get_parameters()->get_vol4d() * NC * (NDIM - 1));
	*splaq /= static_cast<hmc_float>(get_parameters()->get_vol4d() * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	return plaq * 2.0 / static_cast<hmc_float>(get_parameters()->get_vol4d() * NDIM * (NDIM - 1) * NC);
}

hmc_complex Gaugefield_hybrid::polyakov()
{
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int n = 0; n < parameters->get_volspace(); n++) {
		Matrixsu3 prod;
		prod = local_polyakov(get_sgf(), n, parameters);
		hmc_complex tmpcomplex = trace_Matrixsu3(prod);
		complexaccumulate(&res, &tmpcomplex);
	}

	res.re /= static_cast<hmc_float>(NC * parameters->get_volspace());
	res.im /= static_cast<hmc_float>(NC * parameters->get_volspace());
	return res;
}

