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

#ifndef TESTUTILS_H_
#define TESTUTILS_H_

//todo: remove this file eventually

#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"

std::string defaultGpuOption_tmp = "--use_gpu=false";
std::string defaultRec12Option_tmp = "--use_rec12=false";
std::string defaultSourceDirectory_tmp = "../../tests";

static void setArguments(std::string & inputfile_location, std::string & gpu_opt, std::string & rec12_opt, int & num_par, const int param_expect)
{
	logger.info() << "expect command line parameters:";
  logger.info() << "\t<exec_name>\t<source-dir>\t<gpu_usage>\t<rec12_usage>";
	
	switch(num_par){
		case 0:
			logger.fatal() << "Something went terribly wrong! Did not even get executable name! Aborting...";
			exit(-1);
			break;
		case 1:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t" << defaultSourceDirectory_tmp << "\t" << defaultGpuOption_tmp << "\t" << defaultRec12Option_tmp;
			break;
		case 2:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t<source-dir>\t" << defaultGpuOption_tmp << "\t" << defaultRec12Option_tmp;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			break;
		case 3:
			logger.fatal() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.fatal() << "\t<exec_name>\t<source-dir>\t<gpu-usage>\t" << defaultRec12Option_tmp;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			break;

		default:
			if(num_par > param_expect)
			{
    		logger.warn() << "Got " << num_par << " command line parameters, expected only " << param_expect << "! Use only the first " << param_expect << " values";
				num_par = 4;
			}
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
			break;
	}
}

meta::Inputparameters create_parameters(std::string inputfile)
{
	std::string inputfile_location = defaultSourceDirectory_tmp;
	std::string gpu_opt = defaultGpuOption_tmp;
	std::string rec12_opt = defaultRec12Option_tmp;
	int num_par = 0;
	const int param_expect = 4;
  
	num_par = boost::unit_test::framework::master_test_suite().argc;
	setArguments(inputfile_location, gpu_opt, rec12_opt, num_par, param_expect);
	inputfile_location += '/' + inputfile;
	logger.info() << "inputfile used: " << inputfile_location;
	
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

#endif
