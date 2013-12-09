#include "hmcExecutable.h"

hmcExecutable::hmcExecutable(int argc, const char* argv[]) :
	generationExecutable(argc, argv)
{
	initializationTimer.reset();
	meta::print_info_hmc(ownName, parameters);
	writeHmcLogfile();
	setIterationParameters();
	initializationTimer.add();
}

hmcExecutable::~hmcExecutable()
{
  using namespace std;
  logger.info() << "Acceptance rate: " << fixed <<  setprecision(1) << percent(acceptanceRate, parameters.get_hmcsteps()) << "%";
}

inline void hmcExecutable::writeHmcLogfile()
{
	outputToFile.open(filenameForHmcLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_heatbath(ownName, &outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForHmcLogfile);
	}
}

void hmcExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps += parameters.get_hmcsteps();
}

void hmcExecutable::thermalizeAccordingToSpecificAlgorithm()
{
	logger.warn() << "Thermalization is not yet implemented for HMC algorithm";
}

void hmcExecutable::generateAccordingToSpecificAlgorithm()
{
	const double randomNumber = prng->get_double();
	observables = physics::algorithms::perform_hmc_step(gaugefield, iteration, randomNumber, *prng, *system);
	acceptanceRate += observables.accept;
}

void hmcExecutable::performOnlineMeasurements()
{
	if( ( (iteration + 1) % writeFrequency ) == 0 ) {
		std::string gaugeout_name = meta::get_hmc_obs_file_name(parameters, "");
		printHmcObservables(gaugeout_name);
	}
}

void hmcExecutable::printHmcObservables(const std::string& filename)
{
  printHmcObservablesToFile(filename);
  printHmcObservablesToScreen();
}

void hmcExecutable::printHmcObservablesToFile(const std::string& filename)
{
	const hmc_float exp_deltaH = std::exp(observables.deltaH);
	outputToFile.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile << iteration << "\t";
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << observables.plaq << "\t" << observables.tplaq << "\t" << observables.splaq;
	outputToFile << "\t" << observables.poly.re << "\t" << observables.poly.im << "\t" << sqrt(observables.poly.re * observables.poly.re + observables.poly.im * observables.poly.im);
	outputToFile <<  "\t" << observables.deltaH << "\t" << exp_deltaH << "\t" << observables.prob << "\t" << observables.accept;
	//print number of iterations used in inversions with full and force precision
	/**
	 * @todo: The counters should be implemented once the solver class is used!"
	 * until then, only write "0"!
	 */
	int iter0 = 0;
	int iter1 = 0;
	outputToFile << "\t" << iter0 << "\t" << iter1;
	if(parameters.get_use_mp() ) {
		outputToFile << "\t" << iter0 << "\t" << iter1;
	}
	if(meta::get_use_rectangles(parameters) ) {
		outputToFile << "\t" << observables.rectangles;
	}
	outputToFile << std::endl;
	outputToFile.close();
}

void hmcExecutable::printHmcObservablesToScreen()
{
	logger.info() << "\tHMC [OBS]:\t" << iteration << std::setw(8) << std::setfill(' ') << "\t" << std::setprecision(15) 
		      << observables.plaq << "\t" << observables.poly.re << "\t" << observables.poly.im;
}

int main(int argc, const char* argv[])
{
	try {
		hmcExecutable hmcInstance(argc, argv);
		hmcInstance.generateConfigurations();
	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}
