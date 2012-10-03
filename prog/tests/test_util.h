meta::Inputparameters create_parameters(std::string inputfile)
{
  const int param_expect = 4;
  logger.info() << "expect parameters:";
  logger.info() << "\texec_name\tsource-dir\tgpu_usage\trec12_usage";
  std::string inputfile_location, gpu_opt , rec12_opt;
  //get number of parameters                                                                                                                                                                                                                 
  int num_par = boost::unit_test::framework::master_test_suite().argc;
  if(num_par < param_expect){
    logger.fatal() << "Got only " << num_par << " Inputparameters, expected " << param_expect << "! Use inputfile values instead!";
  } else if(num_par > param_expect){
    logger.warn() << "Got " << num_par << " Inputparameters, expected only " << param_expect << "! Use only the first " << param_expect << " values";
  }
  
  switch(num_par){
  case 0:
    logger.fatal() << "Something went terribly wrong! Did not even get executable name! Aborting...";
    exit(-1);
    break;
  case 1:
    logger.fatal() << "Need at least an source-dir! Aborting...";
    exit(-1);
    break;
  case 2:
    //get input file location that has been passed as an argument
    inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
    inputfile_location += inputfile;
    logger.info() << "inputfile used: " << inputfile_location;
    gpu_opt = "";
    rec12_opt = "";
    break;
  case 3:
    //get input file location that has been passed as an argument
    inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
    inputfile_location += inputfile;
    logger.info() << "inputfile used: " << inputfile_location;
    //get use_gpu = true/false that has been passed as an argument
    gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
    logger.info() << "GPU usage: " << gpu_opt;
    rec12_opt = "";
    logger.info() << rec12_opt;
    break;
  default:
    //get input file location that has been passed as an argument
    inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
    inputfile_location += inputfile;
    logger.info() << "inputfile used: " << inputfile_location;
    //get use_gpu = true/false that has been passed as an argument
    gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
    logger.info() << "GPU usage: " << gpu_opt;
    //get use_rec12 = true/false that has been passed as an argument
    rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
    logger.info() << "rec12 usage: " << rec12_opt;
    break;
  }
  const char* _params_cpu[] = {"foo", inputfile_location.c_str(), gpu_opt.c_str() , rec12_opt.c_str()};
  meta::Inputparameters params(num_par, _params_cpu);
  return params;
}

void printKernelInfo(std::string name)
{
  logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}
