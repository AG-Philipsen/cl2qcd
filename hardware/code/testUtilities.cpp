#include "testUtilities.hpp"

#include "../../host_functionality/logger.hpp"
#include <boost/test/unit_test.hpp>

std::string defaultGpuOption = "--use_gpu=false";
std::string defaultRec12Option = "--use_rec12=false";
std::string defaultSourceDirectory = "../../../hardware/code";
std::string defaultInputfilesSubdirectory = "/inputfiles";

static void setArguments(std::string & inputfile_location, std::string & gpu_opt, std::string & rec12_opt, int & num_par, const int param_expect)
{
  logger.info() << "expect command line parameters:";
  logger.info() << "\t<exec_name>\t<source-dir>\t<gpu_usage>\t<rec12_usage>";

	switch(num_par){
		case 0:
			logger.warn() << "Something went terribly wrong! Did not even get executable name! Aborting...";
			exit(-1);
			break;
		case 1:
			logger.warn() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.warn() << "\t<exec_name>\t" << defaultSourceDirectory << "\t" << defaultGpuOption << "\t" << defaultRec12Option;
			break;
		case 2:
			logger.warn() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.warn() << "\t<exec_name>\t<source-dir>\t" << defaultGpuOption << "\t" << defaultRec12Option;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			inputfile_location += defaultInputfilesSubdirectory;
			break;
		case 3:
			logger.warn() << "Got only " << num_par << " command line parameters, expected " << param_expect << "! Use default values instead:";
			logger.warn() << "\t<exec_name>\t<source-dir>\t<gpu-usage>\t" << defaultRec12Option;
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			inputfile_location += defaultInputfilesSubdirectory;
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			break;

		default:
			if(num_par > param_expect)
			{
			  logger.warn() << "Got " << num_par << " command line parameters, expected only " << param_expect << "! Use only the first " << param_expect << " values";
			  num_par = 4;
			}
			
			inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
			inputfile_location += defaultInputfilesSubdirectory;
			gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
			rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
			break;
	}
}

std::unique_ptr<meta::Inputparameters> createParameters(std::string inputfile)
{
	std::string inputfile_location = defaultSourceDirectory + defaultInputfilesSubdirectory;
	std::string gpu_opt = defaultGpuOption;
	std::string rec12_opt = defaultRec12Option;
	int num_par = 0;
	const int param_expect = 4;
  
	num_par = boost::unit_test::framework::master_test_suite().argc;
	setArguments(inputfile_location, gpu_opt, rec12_opt, num_par, param_expect);
	inputfile_location += '/' + inputfile;
	logger.info() << "inputfile used: " << inputfile_location;
	
	const char* _params_cpu[] = {"foo", inputfile_location.c_str(), gpu_opt.c_str() , rec12_opt.c_str(), "--device=0"};
	return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(num_par + 1, _params_cpu));
}

std::unique_ptr<meta::Inputparameters> createParameters(uint numberOfArguments, const char * parameterStringArray[])
{
	const char **newv = (const char**) malloc((numberOfArguments + 3) * sizeof(newv) );
	
	for (uint i = 0; i< numberOfArguments; i++)
	{
		newv[i] = parameterStringArray[i];
	}
	
	std::string inputfile_location = "";
	std::string gpu_opt = defaultGpuOption;
	std::string rec12_opt = defaultRec12Option;
	int num_par = 0;
	const int param_expect = 4;
  
	num_par = boost::unit_test::framework::master_test_suite().argc;
	setArguments(inputfile_location, gpu_opt, rec12_opt, num_par, param_expect);
	
	newv[numberOfArguments] = gpu_opt.c_str();
	newv[numberOfArguments+1] = rec12_opt.c_str();
	newv[numberOfArguments+2] = "--device=0";
	
	for (uint i = 0; i< numberOfArguments + 3; i++)
	{
		logger.info() << newv[i];
	}
	
 	return std::unique_ptr<meta::Inputparameters>(new meta::Inputparameters(3 + numberOfArguments, newv));
}

void printKernelInformation(std::string name)
{
  logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}

