/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

//settings for workaround
std::string use_gpu = "false";
std::string use_rec12 = "false";
std::string std_sourcedir = "../../tests";

meta::Inputparameters create_parameters(std::string inputfile)
{
	std::string inputfile_location, gpu_opt , rec12_opt;
	int num_par = 0;
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
	const char* _params_cpu[] = {"foo", inputfile_location.c_str(), gpu_opt.c_str() , rec12_opt.c_str(), "--device=0"};
  meta::Inputparameters params(num_par + 1, _params_cpu);
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

void testFloatSizeAgainstInputparameters(hmc_float cpu_res, meta::Inputparameters params)
{
	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_SMALL(cpu_res, ref_val);
}
