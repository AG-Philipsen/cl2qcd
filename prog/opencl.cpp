#include "opencl.h"

using namespace std;

hmc_error opencl::init(cl_device_type wanted_device_type){
  cl_int clerr = CL_SUCCESS;

  cout<<"OpenCL being initialized..."<<endl;
  vector<cl::Device> devices;
  vector<cl::Platform> platforms;
  clerr = cl::Platform::get(&platforms);
  cout<<"Number of platforms:    "<<platforms.size()<<endl;
  //LZ: for now, stick to platform 0 without any further checks...
  cout<<"Choose platform number: 0"<<endl;
  cl::Platform platform = platforms[0];
  string info;
  platform.getInfo(CL_PLATFORM_NAME,&info);
  cout<<"\tCL_PLATFORM_NAME:     "<<info<<endl;
  platform.getInfo(CL_PLATFORM_VENDOR,&info);
  cout<<"\tCL_PLATFORM_VENDOR:   "<<info<<endl;
  platform.getInfo(CL_PLATFORM_VERSION,&info);
  cout<<"\tCL_PLATFORM_VERSION:  "<<info<<endl;
  //    platform.getInfo(CL_PLATFORM_PROFILE,&info);
  //    cout<<"\tCL_PLATFORM_PROFILE:  "<<info<<endl;
  cout<<endl;
  platform.getDevices(CL_DEVICE_TYPE_ALL,&devices);
  cout<<"\tAvailable devices: "<<devices.size()<<endl;
  int ndev = 0;
  for(unsigned int m=0; m<devices.size(); m++) {
    string info;
    devices[m].getInfo(CL_DEVICE_NAME,&info);
    cout<<"\t\tCL_DEVICE_NAME:    "<<info<<endl;
    devices[m].getInfo(CL_DEVICE_VENDOR,&info);
    cout<<"\t\tCL_DEVICE_VENDOR:  "<<info<<endl;
    cl_device_type type;
    devices[m].getInfo(CL_DEVICE_TYPE,&type);
    if(type == CL_DEVICE_TYPE_CPU) cout<<"\t\tCL_DEVICE_TYPE:    CPU"<<endl;
    if(type == CL_DEVICE_TYPE_GPU) cout<<"\t\tCL_DEVICE_TYPE:    GPU"<<endl;
    if(type == CL_DEVICE_TYPE_ACCELERATOR) cout<<"\t\tCL_DEVICE_TYPE:    ACCELERATOR"<<endl;
    if(type != CL_DEVICE_TYPE_CPU && type != CL_DEVICE_TYPE_GPU && type != CL_DEVICE_TYPE_ACCELERATOR) {
      cout<<"unexpected CL_DEVICE_TYPE..."<<endl;
      exit(HMC_OCLERROR);
    }
    if(type == wanted_device_type) ndev = m;
    devices[m].getInfo(CL_DEVICE_VERSION,&info);
    cout<<"\t\tCL_DEVICE_VERSION: "<<info<<endl;
  }
  cout<<"Choose device number:   "<<ndev<<endl;
  cout<<"Create context..."<<endl;
  cl_context_properties contextproperties[3] = {CL_CONTEXT_PLATFORM,(cl_context_properties)(platform)(),0};
  const vector<cl::Device> wanted_device = {devices[ndev]};
  const cl::Context context(wanted_device,contextproperties,NULL,NULL,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout<<"Create command queue..."<<endl;
  //eventually, we want to have several queues, one for each device...
  cl::CommandQueue cmdqueue(context,devices[ndev],0,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  //LZ: define sources here if you want to use a vector and several source files...
  //    that might be useful at some point but I haven't been able to make it work so far...
  //  cl::Program::Sources sources;                                        //<SourceVec>
  stringstream sourcecode;
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
    //cout<<"kernel source code:"<<endl; 
    //cout.write(kernelssource,length);
    //cout<<endl;

    //sources.push_back(make_pair(kernelssource,length));                  //<SourceVec>

    kernelsfile.close();
    sourcecode<<kernelssource<<endl;

    delete [] kernelssource;    
  }
  sourcecode<<endl;
  sourcecode<<"//EOF";

  // print complete source code to file
  ofstream kernelsout;
  kernelsout.open("cl_kernelsource.cl");
  if(kernelsout.is_open()){
    kernelsout<<sourcecode.str().c_str()<<endl;
    kernelsout.close();
  } else {
    cout<<"could not open cl_kernelsource.cl"<<endl;
  }

  cout<<"Create program..."<<endl;
  //const cl::Program::Sources csources = sources;   //<SourceVec>
  cl::Program::Sources csources(1,std::make_pair(sourcecode.str().c_str(),sourcecode.str().size()+1));
  cl::Program clprogram(context,csources,&clerr);
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }

  cout<<"Build program..."<<endl;
  stringstream collect_options;
  collect_options<<"-D_INKERNEL_";
#ifdef _RECONSTRUCT_TWELVE_
  collect_options<<" -D_RECONSTRUCT_TWELVE_";
#endif
  collect_options<<" -I"<<SOURCEDIR;
  string buildoptions = collect_options.str();
  cout<<"\tbuild options:";
  cout<<"\t"<<buildoptions<<endl;
  clprogram.build(wanted_device,buildoptions.c_str());
  if(clerr!=CL_SUCCESS) {
    cout<<"... failed, aborting."<<endl;
    exit(HMC_OCLERROR);
  }
  cout<<"Build Log:"<<endl;
  cout<< clprogram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[ndev]).c_str()<<endl;

  return HMC_SUCCESS;
}

hmc_error opencl::finalize_opencl(){
  //this does not work but I can't find anything about Release methods in the C++ API
  //  clReleaseProgram(clprogram);
  //  clReleaseContext(context);
  return HMC_SUCCESS;
}
