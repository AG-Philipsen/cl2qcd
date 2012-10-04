#define BOOST_FAILURE_WORKAROUND 0

//settings for workaround
std::string use_gpu = "false";
std::string use_rec12 = "false";
std::string std_sourcedir = "../../tests";

meta::Inputparameters create_parameters(std::string inputfile)
{
	std::string inputfile_location, gpu_opt , rec12_opt;
	int num_par = 0;
#if BOOST_FAILURE_WORKAROUND
	//CP: this is a workaround in case the boost library is not able to pass arguments
	//	in this case, always use standard args...
	logger.warn() << "Passing of args with boost does not work! Use standard arguments:";
	inputfile_location = std_sourcedir + inputfile;
	gpu_opt = "--use_gpu=" + use_gpu;
	rec12_opt = "--use_rec12=" + use_rec12;
	num_par = 4;
	logger.warn() << inputfile_location << " " << gpu_opt << " " << rec12_opt;
#else
  const int param_expect = 4;
  logger.info() << "expect parameters:";
  logger.info() << "\texec_name\tsource-dir\tgpu_usage\trec12_usage";
  //get number of parameters                                                                                                                                                                                                                 
  num_par = boost::unit_test::framework::master_test_suite().argc;
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
#endif
	const char* _params_cpu[] = {"foo", inputfile_location.c_str(), gpu_opt.c_str() , rec12_opt.c_str()};
  meta::Inputparameters params(num_par, _params_cpu);
  return params;
}

void printKernelInfo(std::string name)
{
  logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}

void testFloatAgainstInputparameters(hmc_float cpu_res, meta::Inputparameters params)
{
	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
}
