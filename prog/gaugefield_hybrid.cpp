#include "gaugefield_hybrid.h"

void Gaugefield_hybrid::init(int numtasks, cl_device_type primary_device_type, inputparameters* input_parameters) {

  logger.trace()<<"Initialize gaugefield";

  //how many tasks (devices with different purpose) do we need:
  set_num_tasks(numtasks);
  //LZ: for now assume that there is only one device per task
  
  //input parameters
  set_parameters(input_parameters);

  //allocate memory for private gaugefield on host and initialize (cold start, read in, future: hot start)
  sgf = new s_gaugefield[1];
  init_gaugefield();

  init_devicetypearray(primary_device_type);
  init_opencl();

  this->init_tasks();

  //this has to be done anyways...
  copy_gaugefield_to_all_tasks();

  return;

}

void Gaugefield_hybrid::init_devicetypearray(cl_device_type primary_device_type){
  //init devicetype array
  if(get_parameters()->get_use_gpu() == false)
    logger.warn()<<"GPU usage turned off in input parameters. Overruled.";

  //LZ: Note that num_task_types and the input parameter num_dev seem not to fit to each other. However, consider following scenario:
  //    num_task_types is given to this class and can potentially control whether certain additional tasks are performed or not
  //    num_dev controls how many devices should be used/are available (so usually CPU and GPU, but possibly also CPU, GPU from different nodes)
  //    here, the different devices can be assigned automatically to the tasks
	//CP: At the moment, get_num_dev is simply not used, altough it may be in the future...
	if(get_parameters()->get_num_dev() != 2)
    logger.warn()<<"Number of devices set to " << get_parameters()->get_num_dev() <<" in input parameters. Overruled: Number of devices must be 2.";

	
  devicetypes = new cl_device_type[get_num_tasks()];
	if(get_num_tasks() == 1){
		devicetypes[0] = primary_device_type;
	}
	else if (get_num_tasks() == 2){
		devicetypes[0] = primary_device_type;
		if(primary_device_type == CL_DEVICE_TYPE_GPU){
			devicetypes[1] = CL_DEVICE_TYPE_CPU;
		} else {
			devicetypes[1] = CL_DEVICE_TYPE_GPU;
		}
	}
	else{
		throw Print_Error_Message("3 or more tasks not yet implemented.");
	}
  
  return;
}

void Gaugefield_hybrid::init_opencl(){
  cl_int clerr = CL_SUCCESS;

  //Initialize OpenCL,
  logger.trace() << "OpenCL being initialized...";

  cl_uint num_platforms;
  //LZ: for now, stick to one platform without any further checks...
  clerr = clGetPlatformIDs(1, &platform, &num_platforms);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetPlatformIDs",__FILE__,__LINE__);


  //Cout Platforminfo
  char info[512];
  clerr = clGetPlatformInfo(platform, CL_PLATFORM_NAME, 512 * sizeof(char), info, NULL);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetPlatformInfo",__FILE__,__LINE__);
  logger.info() << "\tCL_PLATFORM_NAME:     " << info;
  clerr = clGetPlatformInfo(platform, CL_PLATFORM_VENDOR, 512 * sizeof(char), info, NULL);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetPlatformInfo",__FILE__,__LINE__);
  logger.info() << "\tCL_PLATFORM_VENDOR:   " << info;
  clerr = clGetPlatformInfo(platform, CL_PLATFORM_VERSION, 512 * sizeof(char), info, NULL);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetPlatformInfo",__FILE__,__LINE__);
  logger.info() << "\tCL_PLATFORM_VERSION:  " << info;


  cl_uint num_devices_gpu;
  cl_uint num_devices_cpu;
  clerr = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices_gpu);
  clerr = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &num_devices_cpu);
  logger.info() << "Found " << num_devices_gpu << " GPU(s) and " << num_devices_cpu << " CPU(s).";
	
	//Check if number of devices fits progs needs
	if(get_num_tasks() == 1){
		if (devicetypes[0] == CL_DEVICE_TYPE_GPU && num_devices_gpu < 1) 
			throw Print_Error_Message("Application needs one GPU device.");
		else if (devicetypes[0] == CL_DEVICE_TYPE_CPU && num_devices_cpu < 1) 
			throw Print_Error_Message("Application needs one CPU device.");
	}
	else if (get_num_tasks() == 2){
	  if( num_devices_gpu + num_devices_cpu == 1 ){
	    logger.warn() << "You wanted to have two devices, but only one has been found!";
	    set_num_tasks(1);
	    if( num_devices_gpu == 1 ) devicetypes[0] = CL_DEVICE_TYPE_GPU;
	    if( num_devices_cpu == 1 ) devicetypes[0] = CL_DEVICE_TYPE_CPU;
	  }
	  if( num_devices_gpu == 0 ) logger.warn() << "No GPU found.";
	  if( num_devices_cpu == 0 ) logger.warn() << "No CPU found.";
	  if( num_devices_gpu + num_devices_cpu < 1 ) throw Print_Error_Message("Application needs two devices.");
	}

  queue   = new cl_command_queue [get_num_tasks()];
  devices = new cl_device_id     [get_num_tasks()];
  device_double_extension = new string  [get_num_tasks()];
  max_compute_units       = new cl_uint [get_num_tasks()];
  

  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    logger.info()<<"\tInitialize task #"<<ntask<<":";
    clerr = clGetDeviceIDs(platform, get_device_type(ntask), 1, &devices[ntask], NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceIDs",__FILE__,__LINE__);
    
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_NAME, 512 * sizeof(char), info, NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);
    logger.info() << "\t\t\tCL_DEVICE_NAME:    " << info;
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_VENDOR, 512 * sizeof(char), info, NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);
    logger.info() << "\t\t\tCL_DEVICE_VENDOR:  " << info;
    cl_device_type devtype;
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_TYPE, sizeof(cl_device_type), &devtype, NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);
    if(devtype == CL_DEVICE_TYPE_CPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    CPU";
    if(devtype == CL_DEVICE_TYPE_GPU) logger.info() << "\t\t\tCL_DEVICE_TYPE:    GPU";
    if(devtype == CL_DEVICE_TYPE_ACCELERATOR) logger.info() << "\t\t\tCL_DEVICE_TYPE:    ACCELERATOR";
    if(devtype != CL_DEVICE_TYPE_CPU && devtype != CL_DEVICE_TYPE_GPU && devtype != CL_DEVICE_TYPE_ACCELERATOR)
      throw Print_Error_Message("Unexpected CL_DEVICE_TYPE...",__FILE__,__LINE__);
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_VERSION, 512 * sizeof(char), info, NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);
    logger.info() << "\t\t\tCL_DEVICE_VERSION: " << info;
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_EXTENSIONS, 512 * sizeof(char), info, NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);
    logger.info() << "\t\t\tCL_DEVICE_EXTENSIONS: " << info;
    
    if( strstr( info, "cl_amd_fp64" ) != NULL ) device_double_extension[ntask]="AMD";
    if( strstr( info, "cl_khr_fp64" ) != NULL ) device_double_extension[ntask]="KHR";
      
    // figure out the number of "cores"
    clerr = clGetDeviceInfo(devices[ntask], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units[ntask], NULL);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clGetDeviceInfo",__FILE__,__LINE__);;
    logger.info() << "\t\t\tCL_DEVICE_MAX_COMPUTE_UNITS: " << max_compute_units[ntask];

    }

   
  //Initilize context
  logger.trace() << "Create context...";
  context = clCreateContext(0, get_num_tasks(), devices, 0, 0, &clerr);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateContext",__FILE__,__LINE__);

  //Initilize queues
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    logger.trace() << "Create command queue for task #"<<ntask<<"...";

#ifdef _PROFILING_	
    queue[ntask] = clCreateCommandQueue(context, devices[ntask], CL_QUEUE_PROFILING_ENABLE, &clerr);
#else
    queue[ntask] = clCreateCommandQueue(context, devices[ntask], 0, &clerr);
#endif
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateCommandQueue",__FILE__,__LINE__);
  }	


  //context-wide buffers
  logger.trace()<<"Creating gaugefield buffer...";
  clmem_gaugefield = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(s_gaugefield), 0, &clerr);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clCreateBuffer",__FILE__,__LINE__);

  return;
}


void Gaugefield_hybrid::init_tasks(){

  opencl_modules = new Opencl_Module* [get_num_tasks()];
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
    opencl_modules[ntask] = new Opencl_Module[1];
    opencl_modules[ntask]->init(queue[ntask], &clmem_gaugefield, get_parameters(), max_compute_units[ntask], get_double_ext(ntask));
  }

  return;
}


void Gaugefield_hybrid::finalize(){

  this->finalize_opencl();
  this->delete_variables();

  return;
}

void Gaugefield_hybrid::delete_variables(){

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

  return;
}

void Gaugefield_hybrid::finalize_opencl(){

  cl_int clerr = CL_SUCCESS;

  for(int ntask=0; ntask<get_num_tasks(); ntask++) {
    clerr = clFlush(queue[ntask]);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFlush",__FILE__,__LINE__);
    clerr = clFinish(queue[ntask]);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
  }

  //  tasks->clear_kernels();  
  //  tasks->clear_buffers();

  clerr = clReleaseMemObject(clmem_gaugefield);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseMemObject",__FILE__,__LINE__);

  for(int ntask=0; ntask<get_num_tasks(); ntask++) {
    clerr = clReleaseCommandQueue(queue[ntask]);
    if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseCommandQueue",__FILE__,__LINE__);
  }
  
  clerr = clReleaseContext(context);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clReleaseContext",__FILE__,__LINE__);

  return;
}

void Gaugefield_hybrid::init_gaugefield(){

  if((get_parameters())->get_startcondition() == START_FROM_SOURCE) {
    sourcefileparameters parameters_source;
    
    //hmc_gaugefield for filetransfer, initialize here, because otherwise it is not needed
    hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
    //tmp gauge field
    hmc_float * gaugefield_tmp;
    gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float) * NDIM * NC * NC * NTIME * VOLSPACE);
    parameters_source.readsourcefile(&(get_parameters()->sourcefile)[0], get_parameters()->get_prec(), &gaugefield_tmp);
    copy_gaugefield_from_ildg_format(gftmp, gaugefield_tmp, parameters_source.num_entries_source);
    copy_gaugefield_to_s_gaugefield (get_sgf(), gftmp);
    free(gaugefield_tmp);
    free(gftmp);
  }
  if(get_parameters()->get_startcondition() == COLD_START) {
    set_gaugefield_cold(get_sgf());
  }
  if(get_parameters()->get_startcondition() == HOT_START) {
    set_gaugefield_hot(get_sgf());
  }
  
  return;
}

void Gaugefield_hybrid::set_gaugefield_cold(s_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	Matrixsu3 tmp;
	tmp = unit_matrixsu3();
	(*field)[mu][n][t] = tmp;
      }
    }
  }
  return;
}


//Implement this
void Gaugefield_hybrid::set_gaugefield_hot(s_gaugefield * field) {
  throw Print_Error_Message("Hot start not yet implemented.",__FILE__,__LINE__);
  return;
}

void Gaugefield_hybrid::copy_gaugefield_to_all_tasks(){
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    copy_gaugefield_to_task(ntask);
  }
  return;
}

void Gaugefield_hybrid::copy_gaugefield_to_task(int ntask){

  if(ntask < 0 || ntask > get_num_tasks() ) {
    logger.warn()<<"Index out of range, copy_gaugefield_to_device does nothing.";
    return;
  }

  ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(sizeof(s_gaugefield));
  
  copy_to_ocl_format(host_gaugefield, get_sgf());
  
  cl_int clerr = clEnqueueWriteBuffer(queue[ntask], clmem_gaugefield, CL_TRUE, 0, sizeof(s_gaugefield), host_gaugefield, 0, 0, NULL);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clEnqueueWriteBuffer",__FILE__,__LINE__);
  
  free(host_gaugefield);
  return;
}

void Gaugefield_hybrid::synchronize(int ntask_reference){
  if(ntask_reference < 0 || ntask_reference > get_num_tasks() ) {
    logger.warn()<<"Index out of range, synchronize_gaugefield does nothing.";
    return;
  }
  clFinish(queue[ntask_reference]);
  copy_gaugefield_from_task(ntask_reference);

  for(int ntask=0; ntask<get_num_tasks(); ntask++) {
    clFinish(queue[ntask]);
    if(ntask != ntask_reference) copy_gaugefield_to_task(ntask);
  }
  return;
}


void Gaugefield_hybrid::copy_gaugefield_from_task(int ntask){

  if(ntask < 0 || ntask > get_num_tasks() ) {
    logger.warn()<<"Index out of range, copy_gaugefield_from_device does nothing.";
    return;
  }

  ocl_s_gaugefield* host_gaugefield =  (ocl_s_gaugefield*) malloc(sizeof(s_gaugefield));

  cl_int clerr = clEnqueueReadBuffer(queue[ntask], clmem_gaugefield, CL_TRUE, 0, sizeof(s_gaugefield), host_gaugefield, 0, NULL, NULL);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clEnqueueReadBuffer",__FILE__,__LINE__);

  copy_from_ocl_format(get_sgf(), host_gaugefield);

  free(host_gaugefield);

  return;
}

cl_mem* Gaugefield_hybrid::get_clmem_gaugefield(){
  return &clmem_gaugefield;
}

void Gaugefield_hybrid::set_num_tasks (int num){
	num_tasks = num;
	return;
}

int Gaugefield_hybrid::get_num_tasks (){
	return num_tasks;
}

int Gaugefield_hybrid::get_max_compute_units(int ntask){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("rndarray index out of range",__FILE__,__LINE__); 
  return max_compute_units[ntask];
}

string Gaugefield_hybrid::get_double_ext(int ntask){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("rndarray index out of range",__FILE__,__LINE__); 
  return device_double_extension[ntask];
}

inputparameters * Gaugefield_hybrid::get_parameters (){
	return  parameters;
}

void Gaugefield_hybrid::set_parameters (inputparameters * parameters_val){
	parameters = parameters_val;
	return;
}


s_gaugefield * Gaugefield_hybrid::get_sgf (){
    return sgf;
}
	
void Gaugefield_hybrid::set_sgf (s_gaugefield * sgf_val){
	sgf = sgf_val;
	return;
}


cl_device_type Gaugefield_hybrid::get_device_type(int ntask){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range",__FILE__,__LINE__); 
  return devicetypes[ntask];
}

void Gaugefield_hybrid::copy_gaugefield_to_s_gaugefield (s_gaugefield * sgfo, hmc_gaugefield * gf){
  for (int d=0; d <NDIM; d++){
    for (int n=0; n<VOLSPACE; n++){
      for (int t=0; t<NTIME; t++){
#ifdef _RECONSTRUCT_TWELVE_
	  (*sgfo)[d][n][t].e00 = (*gf) [0][d][n][t];
	  (*sgfo)[d][n][t].e01 = (*gf) [2][d][n][t];
	  (*sgfo)[d][n][t].e02 = (*gf) [4][d][n][t];
	  (*sgfo)[d][n][t].e10 = (*gf) [1][d][n][t];
	  (*sgfo)[d][n][t].e11 = (*gf) [3][d][n][t];
	  (*sgfo)[d][n][t].e12 = (*gf) [5][d][n][t];
#else
	  (*sgfo)[d][n][t].e00 = (*gf)[0][0][d][n][t];
	  (*sgfo)[d][n][t].e01 = (*gf) [0][1][d][n][t];
	  (*sgfo)[d][n][t].e02 = (*gf) [0][2][d][n][t];
	  (*sgfo)[d][n][t].e10 = (*gf) [1][0][d][n][t];
	  (*sgfo)[d][n][t].e11 = (*gf) [1][1][d][n][t];
	  (*sgfo)[d][n][t].e12 = (*gf) [1][2][d][n][t];
	  (*sgfo)[d][n][t].e20 = (*gf) [2][0][d][n][t];
	  (*sgfo)[d][n][t].e21 = (*gf) [2][1][d][n][t];
	  (*sgfo)[d][n][t].e22 = (*gf) [2][2][d][n][t];
#endif

      }
    }
  }
  return;
}

void Gaugefield_hybrid::copy_s_gaugefield_to_gaugefield(hmc_gaugefield * gf, s_gaugefield * sgfo){
  for (int d=0; d <NDIM; d++){
    for (int n=0; n<VOLSPACE; n++){
      for (int t=0; t<NTIME; t++){
#ifdef _RECONSTRUCT_TWELVE_
	  (*gf) [0][d][n][t] = (*sgfo)[d][n][t].e00;
	  (*gf) [2][d][n][t] = (*sgfo)[d][n][t].e01;
	  (*gf) [4][d][n][t] = (*sgfo)[d][n][t].e02;
	  (*gf) [1][d][n][t] = (*sgfo)[d][n][t].e10;
	  (*gf) [3][d][n][t] = (*sgfo)[d][n][t].e11;
	  (*gf) [5][d][n][t] = (*sgfo)[d][n][t].e12;
#else
	  (*gf)[0][0][d][n][t] = (*sgfo)[d][n][t].e00;
	  (*gf)[0][1][d][n][t] = (*sgfo)[d][n][t].e01;
	  (*gf)[0][2][d][n][t] = (*sgfo)[d][n][t].e02;
	  (*gf)[1][0][d][n][t] = (*sgfo)[d][n][t].e10;
	  (*gf)[1][1][d][n][t] = (*sgfo)[d][n][t].e11;
	  (*gf)[1][2][d][n][t] = (*sgfo)[d][n][t].e12;
	  (*gf)[2][0][d][n][t] = (*sgfo)[d][n][t].e20;
	  (*gf)[2][1][d][n][t] = (*sgfo)[d][n][t].e21;
	  (*gf)[2][2][d][n][t] = (*sgfo)[d][n][t].e22;
#endif
      }
    }
  }
  return;
}


void Gaugefield_hybrid::save(int number){
	//LZ: generalize the following to larger numbers, if necessary...
	stringstream strnumber;
	strnumber.fill('0');
	strnumber.width(5);
	strnumber << right << number;
	stringstream outfilename;
	outfilename << "conf." << strnumber.str();
	string outputfile = outfilename.str();
	save(outputfile);
	return;
}


void Gaugefield_hybrid::save(string outputfile){
	ildg_gaugefield * gaugefield_buf;
	gaugefield_buf = (ildg_gaugefield*) malloc(sizeof(ildg_gaugefield));
	int gaugefield_buf_size = sizeof(ildg_gaugefield) / sizeof(hmc_float);

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;

	hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	copy_s_gaugefield_to_gaugefield(gftmp, get_sgf());
	copy_gaugefield_to_ildg_format(gaugefield_buf, gftmp);

	hmc_float plaq = plaquette();

	int number = 0;

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NTIME, get_parameters()->get_prec(), number, plaq, get_parameters()->get_beta(), get_parameters()->get_kappa(), get_parameters()->get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	free(gaugefield_buf);
	free(gftmp);

	return;
}


hmc_float Gaugefield_hybrid::plaquette(){
	hmc_float tdummy;
	hmc_float sdummy;
	//LZ: calculation of tdummy and sdummy is unnecessary for this function, chance to speed up a little bit...
	return plaquette(&tdummy, &sdummy);
}

hmc_float Gaugefield_hybrid::plaquette(hmc_float* tplaq, hmc_float* splaq){
	hmc_float plaq = 0;
	*tplaq = 0;
	*splaq = 0;
	
	//CP: new method that is not working right now since elementary matrix-function using structs are missing on the host
	/*
	Matrixsu3 prod;

	for(int t = 0; t < NTIME; t++) {
	  for(int n = 0; n < VOLSPACE; n++) {
	    for(int mu = 0; mu < NDIM; mu++) {
		for(int nu = 0; nu < mu; nu++) {
		  prod = local_plaquette(get_sgf(), pos, t, mu, nu );
		  hmc_float tmpfloat = trace_matrixsu3(prod).re;
		  plaq += tmpfloat;
		  if(mu == 0 || nu == 0) {
		      	tplaq += tmpfloat;
		  } else {
			splaq += tmpfloat;
		  }
		}
	      }
	  }
	}
*/
	//CP: old method, this should be replaced!!
	//LZ: for now it works because I have inserted the copy_to/from routines...
	//LZ: eventually, someone should implement the "structured operations" for the host
	

 	hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	copy_s_gaugefield_to_gaugefield(gftmp, get_sgf());


for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				for(int nu = 0; nu < mu; nu++) {
					hmc_su3matrix prod;
					local_plaquette(gftmp, &prod, n, t, mu, nu );
					hmc_float tmpfloat = trace_su3matrix(&prod).re;
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
	
 free(gftmp);

	*tplaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1));
	*splaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	return plaq * 2.0 / static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) * NC);
}


hmc_complex Gaugefield_hybrid::polyakov(){

 	hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	copy_s_gaugefield_to_gaugefield(gftmp, get_sgf());

	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int n = 0; n < VOLSPACE; n++) {
		hmc_su3matrix prod;
		local_polyakov(gftmp, &prod, n);
		hmc_complex tmpcomplex = trace_su3matrix(&prod);
		complexaccumulate(&res, &tmpcomplex);
	}

	free(gftmp);

	res.re /= static_cast<hmc_float>(NC * VOLSPACE);
	res.im /= static_cast<hmc_float>(NC * VOLSPACE);
	return res;
}


hmc_complex Gaugefield_hybrid::spatial_polyakov(int dir)
{

 	hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	copy_s_gaugefield_to_gaugefield(gftmp, get_sgf());

	//assuming dir=1,2, or 3
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int x1 = 0; x1 < NSPACE; x1++) {
		for(int x2 = 0; x2 < NSPACE; x2++) {
			for(int t = 0; t < NTIME; t++) {
				hmc_su3matrix prod;
				unit_su3matrix(&prod);
				for(int xpol = 0; xpol < NSPACE; xpol++) {
					hmc_su3matrix tmp;
					int coord[NDIM];
					coord[0] = t;
					coord[dir] = xpol;
					int next = (dir % (NDIM - 1)) + 1;
					coord[next] = x1;
					int nnext = (next % (NDIM - 1)) + 1;
					coord[nnext] = x2;
					int pos = get_nspace(coord);
					get_su3matrix(&tmp, gftmp, pos, t, dir);
					accumulate_su3matrix_prod(&prod, &tmp);
				}
				hmc_complex tmpcomplex = trace_su3matrix(&prod);
				complexaccumulate(&res, &tmpcomplex);
			}
		}
	}

	free(gftmp);

	res.re /= static_cast<hmc_float>(NC * NSPACE * NSPACE * NTIME);
	res.im /= static_cast<hmc_float>(NC * NSPACE * NSPACE * NTIME);
	return res;
}

void Gaugefield_hybrid::print_gaugeobservables(int iter)
{
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	hmc_complex pol = polyakov();
	logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
	return;
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
	return;
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range",__FILE__,__LINE__); 
  hmc_float plaq  = 0;
  hmc_float tplaq = 0;
  hmc_float splaq = 0;
  hmc_complex pol;
  cl_mem gf = *get_clmem_gaugefield();
  opencl_modules[ntask]->gaugeobservables(gf, &plaq, &tplaq, &splaq, &pol);
  logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
  return;
}

void Gaugefield_hybrid::print_gaugeobservables_from_task(int iter, int ntask, std::string filename){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range",__FILE__,__LINE__); 
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
  return;
}

#ifdef _PROFILING_
void Gaugefield_hybrid::print_profiling(std::string filename) {
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		//this is initialized with length 1, meaning one assumes one device per task
		opencl_modules[ntask]->print_profiling(filename);
	}
}
#endif
