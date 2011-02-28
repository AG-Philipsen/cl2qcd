#include "opencl.h"

using namespace std;

hmc_error opencl::init(cl_device_type wanted_device_type, usetimer* timer){
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
  clerr = clGetDeviceIDs(platform,wanted_device_type,NULL,NULL,&num_devices);
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
    cout<<"Read kernel source from file: "<<cl_kernels_file[n]<<endl;
    fstream kernelsfile;
    kernelsfile.open(cl_kernels_file[n].c_str());
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

  cout<<"Build Log:"<<endl;
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
    clerr = clEnqueueNDRangeKernel(queue,heatbath_even,1,0,&local_work_size,&global_work_size,0,0,NULL);
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
    clerr = clEnqueueNDRangeKernel(queue,heatbath_odd,1,0,&local_work_size,&global_work_size,0,0,NULL);
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
    clerr = clEnqueueNDRangeKernel(queue,overrelax_even,1,0,&local_work_size,&global_work_size,0,0,NULL);
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
    clerr = clEnqueueNDRangeKernel(queue,overrelax_odd,1,0,&local_work_size,&global_work_size,0,0,NULL);
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
  clerr = clEnqueueNDRangeKernel(queue,plaquette,1,0,&local_work_size,&global_work_size,0,0,NULL);
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
  clerr = clEnqueueReadBuffer(queue,clmem_tplaq,CL_TRUE,0,sizeof(hmc_float),&tplaq[0],0,NULL,NULL);
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
  const size_t local_work_size  = VOL4D/2;
  const size_t global_work_size = local_work_size;
 
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
  hmc_float norm=global_squarenorm(solver_test_in);
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
