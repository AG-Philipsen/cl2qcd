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
    //        cout<<"kernel source code:"<<endl; 
    //        cout.write(kernelssource,length);
    //        cout<<endl;

    kernelsfile.close();
    sourcecode.append(kernelssource,length);

    //    cout<<sourcecode<<endl;

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
  collect_options<<"-D_INKERNEL_ -DNSPACE="<<NSPACE<<" -DNTIME="<<NTIME<<" -DVOLSPACE="<<VOLSPACE;
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
  clmem_plaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_splaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_tplaq = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_float),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clmem_polyakov = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_complex),0,&clerr);
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

hmc_error opencl::run_heatbath(double beta, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  cl_int clerr=CL_SUCCESS;
  (*timer).reset();
  
  hmc_float tmp = (hmc_float) beta;
  
  //cout << "updating even links" << endl;
  clerr = clSetKernelArg(heatbath_even,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_even,1,sizeof(hmc_float),&tmp);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_even,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(heatbath_even,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,heatbath_even,1,0,&local_work_size,&global_work_size,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
  }
  
  //cout << "updating odd links" << endl;
  clerr = clSetKernelArg(heatbath_odd,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_odd,1,sizeof(hmc_float),&tmp);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath_odd,3,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  for(int i = 0; i<NDIM; i++){
    clerr = clSetKernelArg(heatbath_odd,2,sizeof(int),&i);
    if(clerr!=CL_SUCCESS) {
      cout<<"clSetKernelArg failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
    clerr = clEnqueueNDRangeKernel(queue,heatbath_odd,1,0,&local_work_size,&global_work_size,0,0,NULL);
    if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
      exit(HMC_OCLERROR);
    }
  }
  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::run_overrelax(double beta, const size_t local_work_size, const size_t global_work_size, usetimer* timer){
  cl_int clerr=CL_SUCCESS;
  
  (*timer).reset();
  
  hmc_float tmp = (hmc_float) beta;
  
  //cout << "overrelaxing even links" << endl;
  clerr = clSetKernelArg(overrelax_even,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg1 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_even,1,sizeof(hmc_float),&tmp);
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
  
  //cout << "overrelaxing odd links" << endl;
  clerr = clSetKernelArg(overrelax_odd,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(overrelax_odd,1,sizeof(hmc_float),&tmp);
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
  }
  clFinish(queue);
  (*timer).add();
  return HMC_SUCCESS;
}

hmc_error opencl::gaugeobservables(const size_t local_work_size, const size_t global_work_size, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, usetimer* timer1, usetimer * timer2){
  cl_int clerr=CL_SUCCESS;
  
  hmc_float plaq = 0, splaq = 0, tplaq = 0;
  hmc_complex pol = hmc_complex_zero;

  //set device-values to zero for new measurement
  clerr = clEnqueueWriteBuffer(queue,clmem_plaq,CL_TRUE,0,sizeof(hmc_float),&plaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_splaq,CL_TRUE,0,sizeof(hmc_float),&splaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_tplaq,CL_TRUE,0,sizeof(hmc_float),&tplaq,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueWriteBuffer(queue,clmem_polyakov,CL_TRUE,0,sizeof(hmc_complex),&pol,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  //measure plaquette
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
  (*timer1).reset();
  clerr = clEnqueueNDRangeKernel(queue,plaquette,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
      cout<<"clEnqueueNDRangeKernel failed, aborting..." << clerr <<endl;
      exit(HMC_OCLERROR);
  }
  clFinish(queue);
  (*timer1).add();
  
  //measure polyakovloop
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
  (*timer2).reset();
  clerr = clEnqueueNDRangeKernel(queue,polyakov,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clFinish(queue);
  (*timer2).add();
  
  //read out values
  clerr = clEnqueueReadBuffer(queue,clmem_plaq,CL_TRUE,0,sizeof(hmc_float),&plaq,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_splaq,CL_TRUE,0,sizeof(hmc_float),&splaq,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_tplaq,CL_TRUE,0,sizeof(hmc_float),&tplaq,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_polyakov,CL_TRUE,0,sizeof(hmc_complex),&pol,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  //two plaquette-measurements per thread -> add. factor of 1/2

  
  tplaq /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1));
  splaq /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1)*(NDIM-2))/2. ;
  plaq  /= static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC)/2.;
  
  (*plaq_out) = plaq;
  (*splaq_out)= splaq;
  (*tplaq_out)= tplaq;

  
  pol.re /= static_cast<hmc_float>(NC*VOLSPACE);
  pol.im /= static_cast<hmc_float>(NC*VOLSPACE);
  
  (*pol_out).re = pol.re;
  (*pol_out).im = pol.im;
  
  return HMC_SUCCESS;
}

hmc_error opencl::testing(){
  cl_int clerr=CL_SUCCESS;

  cout<<"Create test kernel..."<<endl;
  cl_kernel testkernel = clCreateKernel(clprogram,"test",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  int nsteps = 10;
  hmc_float beta = 4.2;

  const size_t local_work_size  = VOL4D/2;
  const size_t global_work_size = local_work_size;


  hmc_float check=1;
  int size_1 = 3*1000;
  int size_2 = 3000;
  
  clmem_A = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  clmem_B = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_complex)*9, 0, &clerr);
  
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

  hmc_ocl_gaugefield* gfout = (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));
  cl_mem clmem_gfout = clCreateBuffer(context,CL_MEM_READ_WRITE,sizeof(hmc_gaugefield),0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"testing: create clmem_check-buffer failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
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
  clerr = clSetKernelArg(testkernel,4,sizeof(cl_mem),&clmem_gfout);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 4 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,5,sizeof(cl_mem),&clmem_rndarray);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 5 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,6,sizeof(cl_mem),&clmem_random_field_int);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 6 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,7,sizeof(cl_mem),&clmem_random_field_float);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 7 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,8,sizeof(cl_mem),&clmem_random_field_su2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 8 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,9,sizeof(int),&size_1);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 9 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,10,sizeof(int),&size_2);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 10 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,11,sizeof(cl_mem),&clmem_A);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 11 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,12,sizeof(cl_mem),&clmem_B);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg 12 failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  
  clerr = clEnqueueNDRangeKernel(queue,testkernel,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

  /*
  clerr = clEnqueueReadBuffer(queue,clmem_check,CL_TRUE,0,sizeof(hmc_float),&check,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_gfout,CL_TRUE,0,sizeof(hmc_gaugefield),gfout,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  hmc_complex result[9], input[9];
  clerr = clEnqueueReadBuffer(queue,clmem_A,CL_TRUE,0,sizeof(hmc_complex)*9,&result,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... 3 failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueReadBuffer(queue,clmem_B,CL_TRUE,0,sizeof(hmc_complex)*9,&input,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... 3 failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  
  cout << "input test matrix: "<< endl;
  for(int i = 0; i<NC; i++){
    for(int j = 0; j<NC; j++){
      cout << "\t(" << (input[j + NC*i]).re << "," << (input[j + NC*i]).im << ")";
    }cout << endl;}cout << endl;
  
  cout << "result of one overrelaxing step: "<< endl;
  for(int i = 0; i<NC; i++){
    for(int j = 0; j<NC; j++){
      cout << "\t(" << (result[j + NC*i]).re << "," << (result[j + NC*i]).im << ")";
    }cout << endl;}cout << endl;
  */
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
  
  cout<<"\ttest functions: result: "<<check<<endl;
  
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

  if(clReleaseMemObject(clmem_check)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_gfout)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_float)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_int)!=CL_SUCCESS) exit(HMC_OCLERROR);

  free(gfout);
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

  if(clReleaseMemObject(clmem_random_field_int)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_float)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_random_field_su2)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_A)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_B)!=CL_SUCCESS) exit(HMC_OCLERROR);
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
