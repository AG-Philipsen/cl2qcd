#include "opencl.h"

using namespace std;

hmc_error opencl::init(cl_device_type wanted_device_type, usetimer* timer){

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

  (*timer).reset();
  cout<<"OpenCL being initialized..."<<endl;
  cl_uint num_platforms;
  cl_platform_id platform;
  //LZ: for now, stick to one platform without any further checks...
  clerr = clGetPlatformIDs(1,&platform,&num_platforms);
  if(clerr!=CL_SUCCESS) {
    cout<<"clGetPlatformIDs failed..."<<endl;
    exit(HMC_OCLERROR);
  }

  char info[512];
  if(clGetPlatformInfo(platform,CL_PLATFORM_NAME,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\tCL_PLATFORM_NAME:     "<<info<<endl;
  if(clGetPlatformInfo(platform,CL_PLATFORM_VENDOR,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\tCL_PLATFORM_VENDOR:   "<<info<<endl;
  if(clGetPlatformInfo(platform,CL_PLATFORM_VERSION,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\tCL_PLATFORM_VERSION:  "<<info<<endl;
  cout<<endl;

  cl_uint num_devices;
  cl_device_id device;
  clerr = clGetDeviceIDs(platform,wanted_device_type,0,NULL,&num_devices);
  if(num_devices==1) {
    cout<<"\t"<<num_devices<<" device of wanted type has been found."<<endl;
  } else {
    cout<<"\t"<<num_devices<<" devices of wanted type have been found. Choosing device number "<<0<<"."<<endl;
  }
  clerr = clGetDeviceIDs(platform,wanted_device_type,1,&device,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clGetDeviceIDs failed..."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"\tDevice information: "<<endl;
  if(clGetDeviceInfo(device,CL_DEVICE_NAME,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\t\tCL_DEVICE_NAME:    "<<info<<endl;
  if(clGetDeviceInfo(device,CL_DEVICE_VENDOR,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\t\tCL_DEVICE_VENDOR:  "<<info<<endl;
  cl_device_type type;
  if(clGetDeviceInfo(device,CL_DEVICE_TYPE,sizeof(cl_device_type),&type,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  if(type == CL_DEVICE_TYPE_CPU) cout<<"\t\tCL_DEVICE_TYPE:    CPU"<<endl;
  if(type == CL_DEVICE_TYPE_GPU) cout<<"\t\tCL_DEVICE_TYPE:    GPU"<<endl;
  if(type == CL_DEVICE_TYPE_ACCELERATOR) cout<<"\t\tCL_DEVICE_TYPE:    ACCELERATOR"<<endl;
  if(type != CL_DEVICE_TYPE_CPU && type != CL_DEVICE_TYPE_GPU && type != CL_DEVICE_TYPE_ACCELERATOR) {
    cout<<"unexpected CL_DEVICE_TYPE..."<<endl;
    exit(HMC_OCLERROR);
  }
  if(clGetDeviceInfo(device,CL_DEVICE_VERSION,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\t\tCL_DEVICE_VERSION: "<<info<<endl;
  if(clGetDeviceInfo(device,CL_DEVICE_EXTENSIONS,512*sizeof(char),info,NULL) != CL_SUCCESS) exit(HMC_OCLERROR);
  cout<<"\t\tCL_DEVICE_EXTENSIONS: "<<info<<endl;

  cout<<"Create context..."<<endl;
  context = clCreateContext(0,1,&device,0,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Create command queue..."<<endl;
  queue = clCreateCommandQueue(context,device,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  string sourcecode;
  for(unsigned int n=0; n<cl_kernels_file.size(); n++){
	stringstream tmp;
	tmp << SOURCEDIR << '/' << cl_kernels_file[n];
    cout<<"Read kernel source from file: "<<tmp.str()<<endl;

    fstream kernelsfile;
    kernelsfile.open(tmp.str().c_str());
    if(!kernelsfile.is_open()) {
      cout<<"Could not open kernels file. Aborting..."<<endl;
      exit(HMC_FILEERROR);
    }
    
    kernelsfile.seekg(0,ios::end);
    int length = kernelsfile.tellg();
    kernelsfile.seekg(0,ios::beg);
    
    char* kernelssource = new char[length];
    
    kernelsfile.read(kernelssource,length);

    kernelsfile.close();
    sourcecode.append(kernelssource,length);

    delete [] kernelssource;    
  }

  string end = "\n//EOF";
  sourcecode.append(end.c_str(),end.size()+1);

  // print complete source code to file
  ofstream kernelsout;
  kernelsout.open("cl_kernelsource.cl");
  if(kernelsout.is_open()){
    kernelsout<<sourcecode.c_str()<<endl;
    kernelsout.close();
  } else {
    cout<<"could not open cl_kernelsource.cl"<<endl;
  }

  cout<<"Create program..."<<endl;
  size_t sourcesize = sourcecode.size()+1;
  char* source = new char[sourcesize];
  strcpy(source,sourcecode.c_str());
  clprogram = clCreateProgramWithSource(context,1,(const char**)&source,&sourcesize,&clerr);
  delete [] source;
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Build program..."<<endl;
  stringstream collect_options;
  collect_options<<"-D_INKERNEL_ -DNSPACE="<<NSPACE<<" -DNTIME="<<NTIME<<" -DVOLSPACE="<<VOLSPACE <<" -DSPINORSIZE="<<SPINORSIZE <<" -DHALFSPINORSIZE="<<HALFSPINORSIZE <<" -DSPINORFIELDSIZE="<<SPINORFIELDSIZE <<" -DEOPREC_SPINORFIELDSIZE="<<EOPREC_SPINORFIELDSIZE;

#ifdef _RECONSTRUCT_TWELVE_
  collect_options<<" -D_RECONSTRUCT_TWELVE_";
#endif
#ifdef _USEDOUBLEPREC_
  collect_options<<" -D_USEDOUBLEPREC_";
#endif
  collect_options<<" -I"<<SOURCEDIR;
  string buildoptions = collect_options.str();
  cout<<"\tbuild options:";
  cout<<"\t"<<buildoptions<<endl;
  clerr = clBuildProgram(clprogram,1,&device,buildoptions.c_str(),0,0);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, but look at BuildLog and abort then."<<endl;
  }

  cout<<"finished building program" << endl << "Build Log:"<<endl;
  size_t logSize;
  clerr |= clGetProgramBuildInfo(clprogram,device,CL_PROGRAM_BUILD_LOG,0,NULL,&logSize);
  if(logSize!=1) {
    char* log = new char[logSize];
    clerr |= clGetProgramBuildInfo(clprogram,device,CL_PROGRAM_BUILD_LOG,logSize,log,NULL);
    cout<<log<<endl<<endl;
    delete [] log;
  }
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Create buffer for gaugefield..."<<endl;
  clmem_gaugefield = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_gaugefield),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout<<"Create buffer for random numbers..."<<endl;
  clmem_rndarray = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_rndarray),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Create buffer for gaugeobservables..."<<endl;
  clmem_plaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float)*global_work_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_splaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float)*global_work_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_tplaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float)*global_work_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_polyakov = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_complex)*global_work_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Create heatbath kernels..."<<endl;
  heatbath_even = clCreateKernel(clprogram,"heatbath_even",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  heatbath_odd = clCreateKernel(clprogram,"heatbath_odd",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  overrelax_even = clCreateKernel(clprogram,"overrelax_even",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  overrelax_odd = clCreateKernel(clprogram,"overrelax_odd",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout<<"Create gaugeobservables kernels..."<<endl;
  plaquette = clCreateKernel(clprogram,"plaquette",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  polyakov = clCreateKernel(clprogram,"polyakov",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  isinit = 1;

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::copy_gaugefield_to_device(hmc_gaugefield* gaugefield, usetimer* timer){
  cout<<"Copy gaugefield to device..."<<endl;
  (*timer).reset();
  hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

  copy_to_ocl_format(host_gaugefield,gaugefield);

  int clerr = clEnqueueWriteBuffer(queue,clmem_gaugefield,CL_TRUE,0,sizeof(hmc_gaugefield),host_gaugefield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  free(host_gaugefield);

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::copy_rndarray_to_device(hmc_rndarray rndarray, usetimer* timer){
  cout<<"Copy randomarray to device..."<<endl;
  (*timer).reset();

  int clerr = clEnqueueWriteBuffer(queue,clmem_rndarray,CL_TRUE,0,sizeof(hmc_rndarray),rndarray,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::get_gaugefield_from_device(hmc_gaugefield* gaugefield, usetimer* timer){
  cout<<"Get gaugefield from device..."<<endl;
  (*timer).reset();
  hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

  int clerr = clEnqueueReadBuffer(queue,clmem_gaugefield,CL_TRUE,0,sizeof(hmc_gaugefield),host_gaugefield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  copy_from_ocl_format(gaugefield,host_gaugefield);

  free(host_gaugefield);

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer){
  cout<<"Get randomarray from device..."<<endl;
  (*timer).reset();

  int clerr = clEnqueueReadBuffer(queue,clmem_rndarray,CL_TRUE,0,sizeof(hmc_rndarray),rndarray,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::run_heatbath(hmc_float beta, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  cl_int clerr=CL_SUCCESS;
  (*timer).reset();
  
  const size_t * local_work_size_p = (local_work_size == 0) ? 0 : &local_work_size;
  
  clerr = clSetKernelArg(heatbath_even,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg0 at heatbath_even failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_even,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg1 at heatbath_even failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_even,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg3 at heatbath_even failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(heatbath_even,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg2 at heatbath_even failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,heatbath_even,1,0,&global_work_size,local_work_size_p,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clFinish(queue);
  }
  
  clerr = clSetKernelArg(heatbath_odd,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg0 at heatbath_odd failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_odd,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg1 at heatbath_odd failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_odd,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg3 at heatbath_odd failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(heatbath_odd,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg2 at heatbath_odd failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,heatbath_odd,1,0,&global_work_size,local_work_size_p,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clFinish(queue);
  }
  (*timer).add();
  return HMC_SUCCESS;
  
}

hmc_error opencl::run_overrelax(hmc_float beta, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  cl_int clerr=CL_SUCCESS;
  
  (*timer).reset();
  
  const size_t * local_work_size_p = (local_work_size == 0) ? 0 : &local_work_size;

  clerr = clSetKernelArg(overrelax_even,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_even,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_even,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(overrelax_even,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg4 failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,overrelax_even,1,0,&global_work_size,local_work_size_p,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..." << clerr <<endl;
      exit(HMC_OCLERROR);
    }
    clFinish(queue);
  }
  
  clerr = clSetKernelArg(overrelax_odd,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_odd,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_odd,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg7 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(overrelax_odd,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg8 failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,overrelax_odd,1,0,&global_work_size,local_work_size_p,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clFinish(queue);
  }
  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::gaugeobservables(const size_t local_work_size, const size_t global_work_size, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, usetimer* timer1, usetimer * timer2){
  cl_int clerr=CL_SUCCESS;
  
  //measure plaquette
  (*timer1).reset();

  hmc_float * plaq= new hmc_float [global_work_size];
  hmc_float * splaq = new hmc_float[global_work_size];
  hmc_float * tplaq = new hmc_float[global_work_size];
  for(int i = 0; i<(int)global_work_size; i++){
    plaq[i] = 0.;
    splaq[i] = 0.;
    tplaq[i] = 0.;
  }
  // FIXME
  const size_t * local_work_size_p = (local_work_size == 0) ? 0 : &local_work_size;

  //set device-values to zero for new measurement
  clerr = clEnqueueWriteBuffer(queue,clmem_plaq,CL_TRUE,0,sizeof(hmc_float)*global_work_size,&plaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_splaq,CL_TRUE,0,sizeof(hmc_float)*global_work_size,&splaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_tplaq,CL_TRUE,0,sizeof(hmc_float)*global_work_size,&tplaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  clerr = clSetKernelArg(plaquette,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(plaquette,1,sizeof(cl_mem),&clmem_plaq);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(plaquette,2,sizeof(cl_mem),&clmem_tplaq);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(plaquette,3,sizeof(cl_mem),&clmem_splaq);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,plaquette,1,0,&global_work_size,local_work_size_p,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..." << clerr <<endl;
      exit(HMC_OCLERROR);
  }
  clFinish(queue);

  //read out values
  clerr = clEnqueueReadBuffer(queue,clmem_plaq,CL_TRUE,0,sizeof(hmc_float),&plaq[0],0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_splaq,CL_TRUE,0,sizeof(hmc_float),&splaq[0],0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  //two plaquette-measurements per thread -> add. factor of 1/2
  tplaq[0] /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1));
  splaq[0] /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1)*(NDIM-2))/2. ;
  plaq[0]  /= static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC)/2.;
  
  (*plaq_out) = plaq[0];
  (*splaq_out)= splaq[0];
  (*tplaq_out)= tplaq[0];

  (*timer1).add();
  
  //measure polyakovloop
  (*timer2).reset();
  hmc_complex * pol = new hmc_complex [global_work_size];
  for(int i = 0; i<(int) global_work_size; i++){
    pol[i] = hmc_complex_zero;
  }
 
  //set device-values to zero for new measurement
  clerr = clEnqueueWriteBuffer(queue,clmem_polyakov,CL_TRUE,0,sizeof(hmc_complex)*global_work_size,&pol,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  clerr = clSetKernelArg(polyakov,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(polyakov,1,sizeof(cl_mem),&clmem_polyakov);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,polyakov,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);

  //read out values
  clerr = clEnqueueReadBuffer(queue,clmem_polyakov,CL_TRUE,0,sizeof(hmc_complex),&pol[0],0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  pol[0].re /= static_cast<hmc_float>(NC*VOLSPACE);
  pol[0].im /= static_cast<hmc_float>(NC*VOLSPACE);
  
  (*pol_out).re = pol[0].re;
  (*pol_out).im = pol[0].im;

  (*timer2).add();
  
  delete [] plaq;
  delete [] splaq;
  delete [] tplaq;
  delete [] pol;
  
  return HMC_SUCCESS;
}

#ifdef _TESTING_
hmc_error opencl::testing(hmc_gaugefield * gaugefield){
  cl_int clerr=CL_SUCCESS;

  cout<<"Begin testing OpenCL functions"<<endl;
  
  cout<<"Create test kernel..."<<endl;
  cl_kernel testkernel = clCreateKernel(clprogram,"test",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  int nsteps = 10;
  hmc_float beta = 4.2;

  //CP: this is made for a 4^4 lattice, otherwise one will run into problems with the VOL4D/2 definition!!!
  const size_t local_work_size  = 0;
  const size_t global_work_size = VOL4D/2;
  const size_t * local_work_size_p = (local_work_size == 0) ? 0 : &local_work_size;
 
  //CP: random number test
  hmc_float check=1;
  int size_1 = 3*1000;
  int size_2 = 3000;
  
  clmem_random_field_int = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(int)*size_1,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_random_field_float = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(float)*size_2,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_random_field_su2 = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float)*4,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_check = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_check,CL_TRUE,0,sizeof(hmc_float),&check,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  
  //CP: heatbath test (no reconstruct_twelve!!): set input matrices taken from a host_configuration
  hmc_su3matrix host_heatbath_link;
  (host_heatbath_link[0][0]).re = -0.349937;
  (host_heatbath_link[0][0]).im = -0.474231;
  (host_heatbath_link[0][1]).re = 0.301203;
  (host_heatbath_link[0][1]).im = -0.515532;
  (host_heatbath_link[0][2]).re = 0.544024;
  (host_heatbath_link[0][2]).im = -0.013786;
  (host_heatbath_link[1][0]).re = -0.789803;
  (host_heatbath_link[1][0]).im = -0.009681;
  (host_heatbath_link[1][1]).re = 0.047851;
  (host_heatbath_link[1][1]).im = 0.474619;
  (host_heatbath_link[1][2]).re = -0.083667;
  (host_heatbath_link[1][2]).im = 0.376251;  
  (host_heatbath_link[2][0]).re = 0.136194;
  (host_heatbath_link[2][0]).im = 0.101084 ;  
  (host_heatbath_link[2][1]).re = -0.637513;
  (host_heatbath_link[2][1]).im = -0.097608;  
  (host_heatbath_link[2][2]).re = 0.451216;
  (host_heatbath_link[2][2]).im = 0.593032;  

  hmc_staplematrix host_heatbath_staple;
  (host_heatbath_staple[0][0]).re = -0.043457;
  (host_heatbath_staple[0][0]).im = 2.350362;
  (host_heatbath_staple[0][1]).re = -2.827883;
  (host_heatbath_staple[0][1]).im = -0.147794;
  (host_heatbath_staple[0][2]).re = 0.239220;
  (host_heatbath_staple[0][2]).im = -0.353712;
  (host_heatbath_staple[1][0]).re = 2.094735;
  (host_heatbath_staple[1][0]).im = 2.500042;
  (host_heatbath_staple[1][1]).re = 0.274546;
  (host_heatbath_staple[1][1]).im = -2.518978;
  (host_heatbath_staple[1][2]).re = -1.411352;
  (host_heatbath_staple[1][2]).im = -0.568307;  
  (host_heatbath_staple[2][0]).re = 1.217603;
  (host_heatbath_staple[2][0]).im = 0.463087;  
  (host_heatbath_staple[2][1]).re = -0.337425;
  (host_heatbath_staple[2][1]).im = -0.909408;  
  (host_heatbath_staple[2][2]).re = 2.763673;
  (host_heatbath_staple[2][2]).im = -2.460356; 
  
  //set up rnd-array 
  int heatbath_rnd_array_size = 1e4;
  hmc_float * heatbath_rnd_array = (hmc_float*) malloc(sizeof(hmc_float)*heatbath_rnd_array_size);
  
  Random ran(50000);
  for (int i = 0; i<heatbath_rnd_array_size; i++){
    heatbath_rnd_array[i] = ran.doub();    
  }
  
  hmc_su3matrix host_heatbath_out;
  int heatbath_host_cter;
  //CP: there are 3 different possiblities to test the heatbath, with more or less random numbers
//   testing_heatbath_norandommat_no123(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2);
//   testing_heatbath_no123(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2, heatbath_rnd_array, &heatbath_host_cter);
   testing_heatbath(&host_heatbath_link, &host_heatbath_staple, &host_heatbath_out, 4.2, heatbath_rnd_array, &heatbath_host_cter);
   
  clmem_heatbath_test_link_in = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_staple_in = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_link_out = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_rnd_array = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float)*heatbath_rnd_array_size, 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_heatbath_test_cter = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int), 0, &clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  hmc_complex tmp[9], tmp2[9];
  int i, j;
  for(i=0; i<NC; i++){
    for(j=0; j<NC; j++){
      tmp[j*NC+i].re = host_heatbath_link[i][j].re;
      tmp[j*NC+i].im = host_heatbath_link[i][j].im;
      tmp2[j*NC+i].re = host_heatbath_staple[i][j].re;
      tmp2[j*NC+i].im = host_heatbath_staple[i][j].im;     
      
  }}
  
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_link_in,CL_TRUE,0,sizeof(hmc_complex)*9,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_staple_in,CL_TRUE,0,sizeof(hmc_complex)*9,&tmp2,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_heatbath_test_rnd_array,CL_TRUE,0,sizeof(hmc_float)*heatbath_rnd_array_size,heatbath_rnd_array,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
    
  //CP: solver test
  
  hmc_ocl_gaugefield* gfout = (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));
  copy_to_ocl_format(gfout, gaugefield);
  cl_mem clmem_gfout = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_gaugefield),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
  clerr = clEnqueueWriteBuffer(queue,clmem_gfout,CL_TRUE,0,sizeof(hmc_gaugefield),gfout,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
 
  //prepare simple fermionfield
  hmc_spinor_field solver_test_in[SPINORFIELDSIZE];
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int a=0; a<NSPIN; a++) {
	for(int j=0; j<NC; j++) {
	  fill_with_one(solver_test_in, n, t, a, j);
	}
      }
    }
  }
  hmc_float norm = global_squarenorm(solver_test_in);
  norm = sqrt(norm);
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    solver_test_in[n].re /= norm;
    solver_test_in[n].im /=norm;
  }
  
  int test_correlator_size = sizeof(hmc_complex)*NSPACE;
  int test_spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;

  cl_mem clmem_solver_test_spinor_in = clCreateBuffer(context,CL_MEM_READ_WRITE,test_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_solver_test_spinor_out = clCreateBuffer(context,CL_MEM_READ_WRITE,test_spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  cl_mem clmem_solver_test_correlator = clCreateBuffer(context,CL_MEM_READ_WRITE,test_correlator_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

  clerr = clEnqueueWriteBuffer(queue,clmem_solver_test_spinor_in,CL_TRUE,0,test_spinorfield_size,solver_test_in,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }  
  
  //CP: set kernel arguments
  
  clerr = clSetKernelArg(testkernel,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,2,sizeof(int),&nsteps);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,3,sizeof(cl_mem),&clmem_check);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,4,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,5,sizeof(cl_mem),&clmem_random_field_int);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,6,sizeof(cl_mem),&clmem_random_field_float);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,7,sizeof(cl_mem),&clmem_random_field_su2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 7 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,8,sizeof(int),&size_1);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 8 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,9,sizeof(int),&size_2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 9 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,10,sizeof(cl_mem),&clmem_heatbath_test_link_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 10 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,11,sizeof(cl_mem),&clmem_heatbath_test_staple_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 11 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,12,sizeof(cl_mem),&clmem_heatbath_test_link_out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 12 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,13,sizeof(cl_mem),&clmem_heatbath_test_rnd_array);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 13 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,14,sizeof(cl_mem),&clmem_heatbath_test_cter);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 14 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,15,sizeof(cl_mem),&clmem_gfout);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 15 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,16,sizeof(cl_mem),&clmem_solver_test_spinor_in);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 16 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,17,sizeof(cl_mem),&clmem_solver_test_spinor_out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 17 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,18,sizeof(cl_mem),&clmem_solver_test_correlator);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 18 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  cout << "Enqueue kernel.." << endl;
  clerr = clEnqueueNDRangeKernel(queue,testkernel,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
 
  //CP: random number test
  
  int * random_int = (int*) malloc(sizeof(int)*size_1);
  float * random_float = (float*) malloc(sizeof(float)*size_2);
  hmc_float su2mat[4];

  clerr = clEnqueueReadBuffer(queue,clmem_random_field_int,CL_TRUE,0,sizeof(int)*size_1,random_int,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_random_field_float,CL_TRUE,0,sizeof(float)*size_2,random_float,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_random_field_su2,CL_TRUE,0,sizeof(hmc_float)*4,su2mat,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  
  //print result of test from LZ
  cout<<"\ttest functions: result: "<<check<<endl;
  
  //CP: print results from random number test
  cout << "\ttest random 1,2,3 (2 - 1/3*average of "<<size_1/3 << " per array position)" <<endl;
  float rnd_check_int_1 = 0;
  float rnd_check_int_2 = 0;
  float rnd_check_int_3 = 0;
  for(int i=0;i<size_1/3;i++){
// 	  cout << random_int[i*3] << "  " << random_int[i*3 + 1] << "  " << random_int[i*3 + 2] << " " ;
	rnd_check_int_1 += random_int[i*3 + 0];
	rnd_check_int_2 += random_int[i*3 + 1];
	rnd_check_int_3 += random_int[i*3 + 2];
  }
  rnd_check_int_1 = rnd_check_int_1/(size_1/3);
  rnd_check_int_2 = rnd_check_int_2/(size_1/3);
  rnd_check_int_3 = rnd_check_int_3/(size_1/3);
  cout << "\tcheck: " << 2. - rnd_check_int_1 << "  " << 2. - rnd_check_int_2 << "  " << 2. - rnd_check_int_3 <<  endl;
  cout << "\ttest random floats: (0.5 - average of "<<size_2 << ")" << endl;
  float rnd_check_float = 0;
  for(int i=0;i<size_2;i++){
    rnd_check_float += random_float[i];
  }
  
  rnd_check_float/=size_2;
  cout << "\tcheck: " << rnd_check_float-0.5 << endl;
  cout << "\ttest su2mat:";
  cout << "\t" << su2mat[0] << "  " << su2mat[1] << "  " << su2mat[2] << "  " << su2mat[3] << endl;
  cout << "\t1 - det is: " << 1. - (su2mat[0]*su2mat[0]+ su2mat[1]*su2mat[1] + su2mat[2]*su2mat[2] + su2mat[3]*su2mat[3]) << endl;

  
  //CP: heatbath_test
  clerr = clEnqueueReadBuffer(queue,clmem_heatbath_test_link_out,CL_TRUE,0,sizeof(hmc_complex)*9,tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  int heatbath_device_cter;
  clerr = clEnqueueReadBuffer(queue,clmem_heatbath_test_cter,CL_TRUE,0,sizeof(int),&heatbath_device_cter,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout << "\tresult of heatbath test: "<< endl;
  for(i = 0; i<NC; i++){
    for(j = 0; j<NC; j++){
      cout << "\t(" << (host_heatbath_out[i][j]).re - (tmp[i + NC*j]).re << ","  << 
              (host_heatbath_out[i][j]).im  - (tmp[i + NC*j]).im << ")";
    }cout << endl;}cout << endl;
  cout << "\thost code needed " << heatbath_host_cter << " random numbers" << endl;
  cout << "\tdevice code needed " << heatbath_host_cter - heatbath_device_cter << " more random numbers than host code" << endl;
  
  //CP: solver test
  //print result 
  
  hmc_complex solver_test_correlator[NSPACE];
  cout << "\ton the host the correlator is:" << endl;

  hmc_gaugefield * host_gaugefield;
  host_gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
  set_gaugefield_cold(host_gaugefield);
  simple_correlator(host_gaugefield, 0.125, 0.06, 0., 0., 0.,1000);

  clerr = clEnqueueReadBuffer(queue,clmem_solver_test_correlator,CL_TRUE,0,test_correlator_size,&solver_test_correlator,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "\tthe correlator is: " << endl;
  for(int z=0; z<NSPACE; z++) {
    printf("\t%d\t(%e,%e)\n",z, solver_test_correlator[z].re, solver_test_correlator[z].im);
  }
  
  //CP: finish
  if(clReleaseMemObject(clmem_check)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_float)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_int)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_su2)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_link_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_staple_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_link_out)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_heatbath_test_cter)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_gfout)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_spinor_in)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_spinor_out)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_solver_test_correlator)!=CL_SUCCESS) exit(HMC_OCLERROR);

  free(gfout);
  
  printf("opencl testing finished\n");
  return HMC_SUCCESS;
}
#endif //_TESTING_

hmc_error opencl::finalize(){
  if(isinit==1) {
  if(clFlush(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clFinish(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseKernel(heatbath_even)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(heatbath_odd)!=CL_SUCCESS) exit(HMC_OCLERROR);
  
  if(clReleaseKernel(overrelax_even)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(overrelax_odd)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseKernel(plaquette)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseKernel(polyakov)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseProgram(clprogram)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseMemObject(clmem_gaugefield)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_rndarray)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseMemObject(clmem_plaq)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_tplaq)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_splaq)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_polyakov)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseCommandQueue(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseContext(context)!=CL_SUCCESS) exit(HMC_OCLERROR);

  isinit = 0;
  }
  return HMC_SUCCESS;
}

#ifdef _FERMIONS_
hmc_error opencl::init_solver_variables(inputparameters* parameters, const size_t local_work_size, const size_t global_work_size, usetimer * timer){
	
	(*timer).reset();
		
	cout << "init solver variables..." << endl;
	int clerr = CL_SUCCESS; 
	//!!CP: ?????
	int num_groups;
	if(local_work_size <= global_work_size) num_groups = global_work_size/local_work_size;
	else num_groups = 1;
	cout<<"num_groups: " << num_groups << endl;
	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int complex_size = sizeof(hmc_complex);
	int float_size = sizeof(hmc_float);
	int local_buf_size = complex_size; //!!sizeof(complex)*local_work_size;
	int global_buf_size = complex_size; //!!sizeof(complex)*( num_groups
	int local_buf_size_float = float_size*local_work_size; //!!sizeof(float)*local_work_size;
	int global_buf_size_float = float_size*num_groups; //!!sizeof(float)*( num_groups
	hmc_complex one = hmc_complex_one;
	hmc_complex minusone = hmc_complex_minusone;
	hmc_complex kappa_complex = {(*parameters).get_kappa(), 0.};
	hmc_float tmp;
	
	cout << "\tinit spinorfields..." << endl;
  clmem_inout = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_inout failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_source = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_source failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rn = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rn failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rhat = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rhat failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_v = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_p = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_p failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_s = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_s failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_t = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_t failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_aux = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp = clCreateBuffer(context,CL_MEM_READ_WRITE,spinorfield_size,0,&clerr);;
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_v failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "\tinit complex numbers..." << endl;
	clmem_rho = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rho failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_rho_next = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_rho_next failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_alpha = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_alpha failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_omega = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_omega failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_beta = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_beta failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp1 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_tmp2 = clCreateBuffer(context,CL_MEM_READ_WRITE,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_tmp2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_one = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_one failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_minusone = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_minusone failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_kappa_cmplx = clCreateBuffer(context,CL_MEM_READ_ONLY,complex_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_kappa_cplx failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
	clmem_scalar_product_buf_loc = clCreateBuffer(context,CL_MEM_READ_WRITE,local_buf_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_scalar_product_buf_loc failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_scalar_product_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_scalar_product_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  

  cout << "\tinit float numbers..." << endl;
	clmem_kappa = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_kappa failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
	clmem_theta_fermion = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_theta_fermion failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_mu = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_mu failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_chem_pot_re = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_chem_pot_re failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_chem_pot_im = clCreateBuffer(context,CL_MEM_READ_ONLY,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_chem_pot_im failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_resid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_resid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
	clmem_trueresid = clCreateBuffer(context,CL_MEM_READ_WRITE,float_size,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_trueresid failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  } 
  clmem_global_squarenorm_buf_loc = clCreateBuffer(context,CL_MEM_READ_WRITE,local_buf_size_float,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_global_squarenorm_buf_loc failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clmem_global_squarenorm_buf_glob = clCreateBuffer(context,CL_MEM_READ_WRITE,global_buf_size_float,0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"creating clmem_global_squarenorm_buf_glob failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  
  cout << "write values to buffers..." << endl;
	tmp = (*parameters).get_kappa();
  clerr = clEnqueueWriteBuffer(queue,clmem_kappa,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_kappa failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_theta_fermion();
  clerr = clEnqueueWriteBuffer(queue,clmem_theta_fermion,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_theta_fermion failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_mu();
  clerr = clEnqueueWriteBuffer(queue,clmem_mu,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_mu failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_chem_pot_re();
  clerr = clEnqueueWriteBuffer(queue,clmem_chem_pot_re,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_chem_pot_re failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  tmp = (*parameters).get_chem_pot_im();
  clerr = clEnqueueWriteBuffer(queue,clmem_chem_pot_im,CL_TRUE,0,float_size,&tmp,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_chem_pot_im failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

	clerr = clEnqueueWriteBuffer(queue,clmem_one,CL_TRUE,0,complex_size,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_one failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  } 
  clerr = clEnqueueWriteBuffer(queue,clmem_minusone,CL_TRUE,0,complex_size,&minusone,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_minusone failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_kappa_cmplx,CL_TRUE,0,complex_size,&kappa_complex,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... writing clmem_kappa_cplx failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "Init fermion kernels..." << endl;
	
	M_diag = clCreateKernel(clprogram,"M_diag",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating M_diag kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	dslash = clCreateKernel(clprogram,"dslash",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating dslash kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxpy = clCreateKernel(clprogram,"saxpy",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxpy kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	saxsbypz = clCreateKernel(clprogram,"saxsbypz",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating saxsbypz kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	scalar_product = clCreateKernel(clprogram,"scalar_product",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating scalar_product kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	scalar_product_reduction = clCreateKernel(clprogram,"scalar_product_reduction",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating scalar_product_reduction kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	set_zero_spinorfield = clCreateKernel(clprogram,"set_zero_spinorfield",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating set_zero_spinorfield kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	global_squarenorm = clCreateKernel(clprogram,"global_squarenorm",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating global_squarenorm kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	global_squarenorm_reduction = clCreateKernel(clprogram,"global_squarenorm_reduction",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating global_squarenorm_reduction kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	ratio = clCreateKernel(clprogram,"ratio",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating ratio kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
	product = clCreateKernel(clprogram,"product",&clerr);
	if(clerr!=CL_SUCCESS) {
		cout<<"...creating product kernel failed, aborting."<<endl;
		exit(HMC_OCLERROR);
	}
  
  (*timer).add();
  return HMC_SUCCESS;
}


hmc_error opencl::copy_spinorfield_to_device(hmc_spinor_field* host_spinorfield,  usetimer* timer){
	cout<<"Copy spinorfield to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_source_to_device(hmc_spinor_field* host_source,  usetimer* timer){
	cout<<"Copy source to device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	int clerr = clEnqueueWriteBuffer(queue,clmem_source,CL_TRUE,0,spinorfield_size,host_source,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::get_spinorfield_from_device(hmc_spinor_field* host_spinorfield, usetimer* timer){
//   cout<<"Get spinorfield from device..."<<endl;
  (*timer).reset();

	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
  int clerr = clEnqueueReadBuffer(queue,clmem_inout,CL_TRUE,0,spinorfield_size,host_spinorfield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }

  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int spinorfield_size = sizeof(hmc_complex)*SPINORFIELDSIZE;
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, spinorfield_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	hmc_float tmp;
	clerr = clEnqueueReadBuffer(queue,in,CL_TRUE,0,sizeof(hmc_float),&tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	hmc_complex tmp;
	clerr = clEnqueueReadBuffer(queue,in,CL_TRUE,0,sizeof(hmc_complex),&tmp,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    cout <<"errorcode :" << clerr << endl;
    exit(HMC_OCLERROR);
  }
	(*out) = tmp;
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::M_device(cl_mem in, cl_mem out, cl_mem tmp, cl_mem gaugefield, cl_mem kappa, cl_mem mu, cl_mem theta, cl_mem kappa_complex, cl_mem chem_pot_re, cl_mem chem_pot_im, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(M_diag,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,2,sizeof(cl_mem),&kappa);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(M_diag,3,sizeof(cl_mem),&mu);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,M_diag,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue M_diag kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
	clerr = clSetKernelArg(dslash,0,sizeof(cl_mem),&in); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,1,sizeof(cl_mem),&tmp);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,2,sizeof(cl_mem),&gaugefield);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,3,sizeof(cl_mem),&theta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,4,sizeof(cl_mem),&chem_pot_re);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(dslash,5,sizeof(cl_mem),&chem_pot_im);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,dslash,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue dslash kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	//!! perhaps this can go into an extra calling of saxpy
	clerr = clSetKernelArg(saxpy,0,sizeof(cl_mem),&tmp); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,2,sizeof(cl_mem),&clmem_kappa_cmplx);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,3,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,saxpy,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxpy kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxpy,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,2,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxpy,3,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,saxpy,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxpy kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(saxsbypz,0,sizeof(cl_mem),&x); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,1,sizeof(cl_mem),&y);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,2,sizeof(cl_mem),&z);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,3,sizeof(cl_mem),&alpha);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,4,sizeof(cl_mem),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(saxsbypz,5,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,saxsbypz,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue saxsbypz kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(scalar_product,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(scalar_product,2,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //!!CP: this must be in local memory!!!!
	clerr = clSetKernelArg(scalar_product,3,sizeof(hmc_complex)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 3 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,scalar_product,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(scalar_product_reduction,0,sizeof(cl_mem),&clmem_scalar_product_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(scalar_product_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,scalar_product_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(ratio,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(ratio,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(ratio,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //!!CP:this needs only one kernel!!
	size_t one = 1;
  clerr = clEnqueueNDRangeKernel(queue,ratio,1,0,&one,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue ratio kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
  (*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer){
	(*timer).reset();
	
	int clerr =CL_SUCCESS;
	clerr = clSetKernelArg(product,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(product,1,sizeof(cl_mem),&b);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(product,2,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //!!CP:this needs only one kernel!!
	size_t one = 1;
  clerr = clEnqueueNDRangeKernel(queue,product,1,0,&one,&one,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue product kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	
	int clerr = CL_SUCCESS;
	clerr = clSetKernelArg(global_squarenorm,0,sizeof(cl_mem),&a);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm,1,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  //!!CP: this must be in local memory!!!!
	clerr = clSetKernelArg(global_squarenorm,2,sizeof(hmc_float)*local_work_size,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 2 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }  
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue global_squarenorm kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  clerr = clSetKernelArg(global_squarenorm_reduction,0,sizeof(cl_mem),&clmem_global_squarenorm_buf_glob);      
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clSetKernelArg(global_squarenorm_reduction,1,sizeof(cl_mem),&out);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,global_squarenorm_reduction,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue scalar_product_reduction kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);

	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::set_zero_spinorfield_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	
	clerr = clSetKernelArg(set_zero_spinorfield,0,sizeof(cl_mem),&x);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 0 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
	clerr = clEnqueueNDRangeKernel(queue,set_zero_spinorfield,1,0,&global_work_size,&local_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"enqueue set_zero_spinorfield kernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
	
	(*timer).add();
	return HMC_SUCCESS;
}

hmc_error opencl::copy_complex_device(cl_mem in, cl_mem out, usetimer* timer){
	(*timer).reset();
	int clerr = CL_SUCCESS;
	int complex_size = sizeof(hmc_complex);
	
	clerr = clEnqueueCopyBuffer(queue,in, out, 0,0, complex_size ,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
	}	
	(*timer).add();
	return HMC_SUCCESS;
}
	
hmc_error opencl::bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, const size_t local_work_size, const size_t global_work_size, int cgmax){

	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	hmc_float trueresid;

	//!!CP: I think most of the args in M_device can be saved since kappa... are known global

	for(int iter=0; iter<cgmax; iter++){
		if(iter%iter_refresh==0) {
			set_zero_spinorfield_device(clmem_v, localsize, globalsize, latimer); 
			set_zero_spinorfield_device(clmem_p, localsize, globalsize, latimer);
			//CP: debugging
			set_zero_spinorfield_device(clmem_tmp, localsize, globalsize, latimer);
			
			M_device(clmem_inout, clmem_rn, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion,  clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, localsize, globalsize, Mtimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
			copy_spinor_device(clmem_rn, clmem_rhat, singletimer);

			
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			copy_complex_device(clmem_one, clmem_omega, singletimer);
			copy_complex_device(clmem_one, clmem_rho, singletimer);
			
			//!!CP: calc initial residuum, this is not needed for the algorithm!!
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_resid, &resid, copytimer);
			cout << "initial residuum is: " << resid << endl;
		}

		set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
		copy_complex_device(clmem_rho_next, clmem_rho, singletimer);
		set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
		set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);

		set_complex_to_product_device(clmem_beta, clmem_omega, clmem_tmp1, singletimer);
		set_complex_to_product_device(clmem_minusone, clmem_tmp1, clmem_tmp2, singletimer);
		saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p, local_work_size, global_work_size, latimer);

		M_device(clmem_p,clmem_v, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion, clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, local_work_size, global_work_size, Mtimer);

		set_complex_to_scalar_product_device(clmem_rhat, clmem_v, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device (clmem_rho, clmem_tmp1, clmem_alpha, singletimer);
		
		saxpy_device(clmem_v, clmem_rn, clmem_alpha, clmem_s, local_work_size, global_work_size, latimer);
		
		M_device(clmem_s, clmem_t, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion, clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, local_work_size, global_work_size, Mtimer);

		set_complex_to_scalar_product_device(clmem_t,clmem_s, clmem_tmp1, local_work_size, global_work_size, scalarprodtimer);
		//!!CP: can this also be global_squarenorm??
		set_complex_to_scalar_product_device(clmem_t,clmem_t, clmem_tmp2, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_tmp1, clmem_tmp2, clmem_omega, singletimer);

		saxpy_device(clmem_t, clmem_s, clmem_omega, clmem_rn, local_work_size, global_work_size, latimer);

		saxsbypz_device(clmem_p, clmem_s, clmem_inout, clmem_alpha, clmem_omega, clmem_inout, local_work_size, global_work_size, latimer); 
    
		set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid<epssquare) {	
			M_device(clmem_inout,clmem_aux,clmem_tmp, clmem_gaugefield,clmem_kappa,clmem_mu, clmem_theta_fermion, clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, local_work_size, global_work_size, Mtimer);
			saxpy_device(clmem_aux, clmem_source, clmem_one, clmem_aux, local_work_size, global_work_size, latimer); 
			set_float_to_global_squarenorm_device(clmem_aux, clmem_trueresid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_trueresid, &trueresid, copytimer);
			cout << "residuum:\t" << resid << "\ttrueresiduum:\t" << trueresid << endl;
			if(trueresid<epssquare)
				return HMC_SUCCESS;
		}
		else{
			cout << "residuum:\t" << resid << endl;
		}

	}

	return HMC_SUCCESS;
}

hmc_error opencl::cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, const size_t local_work_size, const size_t global_work_size, int cgmax){
	//!!CP: here one has to be careful if local_work_size is a null-pointer
	size_t globalsize = global_work_size;
	size_t localsize = local_work_size;
	//CP: these have to be on the host
	hmc_float resid;
	int iter;
  for(iter = 0; iter < cgmax; iter ++){  
    if(iter%iter_refresh == 0){
			M_device(clmem_inout, clmem_rn, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion,  clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, localsize, globalsize, Mtimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
      copy_spinor_device(clmem_rn, clmem_p, copytimer);
    }
    //alpha = (rn, rn)/(pn, Apn) --> alpha = omega/rho
		set_complex_to_scalar_product_device(clmem_p, clmem_p, clmem_omega, local_work_size, global_work_size, scalarprodtimer);
		//A pn --> v
		M_device(clmem_p,clmem_v, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion, clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, local_work_size, global_work_size, Mtimer);
		set_complex_to_scalar_product_device(clmem_p, clmem_v, clmem_rho, local_work_size, global_work_size, scalarprodtimer);
		set_complex_to_ratio_device(clmem_omega, clmem_rho, clmem_alpha, singletimer);
		
		set_complex_to_product_device(clmem_alpha, clmem_minusone, clmem_tmp1, singletimer);

		saxpy_device(clmem_inout, clmem_p, clmem_tmp1, clmem_inout, localsize, globalsize, latimer);
		//rn+1 -> rhat
		saxpy_device(clmem_rn, clmem_v, clmem_alpha, clmem_rhat, localsize, globalsize, latimer);
		
		set_float_to_global_squarenorm_device(clmem_rhat, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
		copy_float_from_device(clmem_resid, &resid, copytimer);

		if(resid<epssquare) {
			return HMC_SUCCESS;
		}
		else{
			//beta = (rn+1, rn+1)/(rn, rn) --> alpha = rho_next/omega
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
			set_complex_to_ratio_device(clmem_rho_next, clmem_omega, clmem_beta, singletimer);
			
			//pn+1 = rn+1 + beta*pn 
			set_complex_to_product_device(clmem_beta, clmem_minusone, clmem_tmp2, singletimer);
			saxpy_device(clmem_p, clmem_rhat, clmem_tmp2, clmem_p, localsize, globalsize, latimer);
		}
	}
		return HMC_SUCCESS;
}
	

//!!CP: this is here because the kernel global_squarenorm has the same name as the function in operations_spinor
// hmc_float global_squarenorm(hmc_spinor_field *field) {
hmc_float inline global_squarenorm_host(hmc_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<NTIME; t++) {
    for (int n=0; n<VOLSPACE; n++) {
      sum += local_squarenorm(field,n,t);
    }
  }
  return sum;
}
	
hmc_error opencl::testing_spinor(size_t local_size, size_t global_size){

	usetimer noop;
	
	cout << "testing spinor kernels..." << endl;
	 //set up spinor field for testing                                                                                                                          
	hmc_spinor_field test[SPINORFIELDSIZE];
  for (int i = 0; i<SPINORFIELDSIZE; i++){
    test[i].re = 1.;
    test[i].im = 0.;
  }
  hmc_float norm=global_squarenorm_host(test);
  norm = sqrt(norm);
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    test[n].re /= norm;
    test[n].im /=norm;
  }
  cout << "\tspinorfield norm: "<< global_squarenorm_host(test) << endl;
  copy_spinorfield_to_device(test, &noop);
  get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after copying to device: "<< global_squarenorm_host(test) << endl;
  copy_spinor_device(clmem_inout, clmem_v, &noop);
	set_zero_spinorfield_device(clmem_inout, local_size, global_size, &noop);
	get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after set_zero_spinorfield: "<< global_squarenorm_host(test) << endl;
	copy_spinor_device(clmem_v, clmem_inout, &noop);
	get_spinorfield_from_device(test, &noop);
  cout << "\tspinorfield norm after copying spinorfield on device: "<< global_squarenorm_host(test) << endl;
	copy_float_from_device(clmem_resid, &norm, &noop);
	cout<< "\tglobal squarenorm on the device before calculation is: " << norm << endl;
	set_float_to_global_squarenorm_device(clmem_inout, clmem_resid, local_size, global_size, &noop);
	copy_float_from_device(clmem_resid, &norm, &noop);
	cout<< "\tglobal squarenorm on the device is: " << norm << endl;
	hmc_complex tester;
	copy_complex_from_device(clmem_tmp1, &tester, &noop);
	cout<< "\ttmp1 before copying is: " << tester.re << "," << tester.im << endl;	
	copy_complex_device(clmem_minusone, clmem_tmp1, &noop);
	copy_complex_from_device(clmem_tmp1, &tester, &noop);
	cout<< "\ttmp1 after copying minusone on it is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_ratio_device(clmem_tmp1, clmem_one, clmem_tmp2, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\ttmp2 after copying is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_product_device(clmem_tmp1, clmem_minusone, clmem_tmp2, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\ttmp2 after copying is: " << tester.re << "," << tester.im << endl;	

	set_complex_to_scalar_product_device(clmem_inout, clmem_v, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(inout,inout) is: " << tester.re << "," << tester.im << endl;	
	
	saxpy_device(clmem_inout, clmem_v, clmem_minusone, clmem_rhat, local_size, global_size, &noop);
  set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t|(inout + inout)| is: " << tester.re << "," << tester.im << endl;
	
	saxsbypz_device(clmem_inout, clmem_v, clmem_inout, clmem_minusone, clmem_one, clmem_rhat, local_size, global_size, &noop);
  set_complex_to_scalar_product_device(clmem_rhat, clmem_rhat, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t|(inout + inout - inout)| is: " << tester.re << "," << tester.im << endl;
	
	M_device(clmem_inout, clmem_v, clmem_tmp, clmem_gaugefield, clmem_kappa, clmem_mu, clmem_theta_fermion, clmem_kappa_cmplx, clmem_chem_pot_re, clmem_chem_pot_im, local_size, global_size, &noop);
	
	set_complex_to_scalar_product_device(clmem_v, clmem_v, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(Mv, Mv) is: " << tester.re << "," << tester.im << endl;	
	
	hmc_spinor_field b[SPINORFIELDSIZE];
	//CP: a gaugefield is not needed here
	hmc_gaugefield dummy;
	create_point_source(b,1,0,0,0.15,4.,&dummy);
	copy_source_to_device(b, &noop);
	set_complex_to_scalar_product_device(clmem_source, clmem_source, clmem_tmp2, local_size, global_size, &noop);
	copy_complex_from_device(clmem_tmp2, &tester, &noop);
	cout<< "\t(source, source) is: " << tester.re << "," << tester.im << endl;	
	
	cout << "perform bicgstab..." << endl;
	bicgstab_device(&noop, &noop, &noop, &noop, &noop,local_size, global_size, 100);
	
	cout << "...testing fermion kernels done" << endl;
	
	return HMC_SUCCESS;
}
	
	
#endif //_FERMIONS_
