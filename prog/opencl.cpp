#include "opencl.h"

using namespace std;

hmc_error opencl::init(cl_device_type wanted_device_type){
  cl_int clerr = CL_SUCCESS;

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
  collect_options<<" -D_USEDOUBLEPREC";
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

  cout<<"Create heatbath kernel..."<<endl;
  heatbath = clCreateKernel(clprogram,"heatbath",&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  isinit = 1;

  return HMC_SUCCESS;
}

hmc_error opencl::copy_gaugefield_to_device(hmc_gaugefield* gaugefield){
  cout<<"Copy gaugefield to device..."<<endl;

  hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

  copy_to_ocl_format(host_gaugefield,gaugefield);

  int clerr = clEnqueueWriteBuffer(queue,clmem_gaugefield,CL_TRUE,0,sizeof(hmc_gaugefield),host_gaugefield,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
		     
  free(host_gaugefield);

  return HMC_SUCCESS;
}

hmc_error opencl::get_gaugefield_from_device(hmc_gaugefield* gaugefield){
  cout<<"Get gaugefield from device..."<<endl;

  hmc_ocl_gaugefield* host_gaugefield =  (hmc_ocl_gaugefield*) malloc(sizeof(hmc_gaugefield));

  int clerr = clEnqueueReadBuffer(queue,clmem_gaugefield,CL_TRUE,0,sizeof(hmc_gaugefield),host_gaugefield,0,NULL,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  copy_from_ocl_format(gaugefield,host_gaugefield);

  free(host_gaugefield);

  return HMC_SUCCESS;
}

hmc_error opencl::run_heatbath(int nsteps, double beta){
  cout<<"perform "<<nsteps<<" heatbath steps on OpenCL device..."<<endl;
  const size_t local_work_size  = VOLSPACE/2;
  const size_t global_work_size = local_work_size;
  cl_int clerr=CL_SUCCESS;
  clerr = clSetKernelArg(heatbath,0,sizeof(cl_mem),&clmem_gaugefield); 
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(heatbath,2,sizeof(int),&nsteps);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,heatbath,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

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

  const size_t local_work_size  = VOLSPACE/2;
  const size_t global_work_size = local_work_size;


  hmc_float check=1;
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
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,1,sizeof(hmc_float),&beta);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,2,sizeof(int),&nsteps);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,3,sizeof(cl_mem),&clmem_check);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clSetKernelArg(testkernel,4,sizeof(cl_mem),&clmem_gfout);
  if(clerr!=CL_SUCCESS) {
    cout<<"clSetKernelArg failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }
  clerr = clEnqueueNDRangeKernel(queue,testkernel,1,0,&local_work_size,&global_work_size,0,0,NULL);
  if(clerr!=CL_SUCCESS) {
    cout<<"clEnqueueNDRangeKernel failed, aborting..."<<endl;
    exit(HMC_OCLERROR);
  }

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

  cout<<"check result: "<<check<<endl;
  if(clReleaseMemObject(clmem_check)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseMemObject(clmem_gfout)!=CL_SUCCESS) exit(HMC_OCLERROR);
  
  free(gfout);

  return HMC_SUCCESS;
}

hmc_error opencl::finalize(){
  if(isinit==1) {
  if(clFlush(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clFinish(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseKernel(heatbath)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseProgram(clprogram)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseMemObject(clmem_gaugefield)!=CL_SUCCESS) exit(HMC_OCLERROR);

  if(clReleaseCommandQueue(queue)!=CL_SUCCESS) exit(HMC_OCLERROR);
  if(clReleaseContext(context)!=CL_SUCCESS) exit(HMC_OCLERROR);

  isinit = 0;
  }
  return HMC_SUCCESS;
}
