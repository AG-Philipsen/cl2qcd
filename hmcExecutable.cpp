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
		print_hmcobservables(gaugeout_name, parameters);
	} else if(parameters.get_print_to_screen() ) {
		print_hmcobservables();
	}
}

void hmcExecutable::print_hmcobservables(const std::string& filename, const meta::Inputparameters& params)
{
	const hmc_float exp_deltaH = std::exp(observables.deltaH);
	std::fstream hmcout(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) throw File_Exception(filename);
	hmcout << iteration << "\t";
	hmcout.width(8);
	hmcout.precision(15);
	//print plaquette (plaq, tplaq, splaq)
	hmcout << observables.plaq << "\t" << observables.tplaq << "\t" << observables.splaq;
	//print polyakov loop (re, im, abs)
	hmcout << "\t" << observables.poly.re << "\t" << observables.poly.im << "\t" << sqrt(observables.poly.re * observables.poly.re + observables.poly.im * observables.poly.im);
	//print deltaH, exp(deltaH), acceptance-propability, accept (yes or no)
	hmcout <<  "\t" << observables.deltaH << "\t" << exp_deltaH << "\t" << observables.prob << "\t" << observables.accept;
	//print number of iterations used in inversions with full and force precision
	/**
	 * @todo: The counters should be implemented once the solver class is used!"
	 * until then, only write "0"!
	 */
	int iter0 = 0;
	int iter1 = 0;
	hmcout << "\t" << iter0 << "\t" << iter1;
	if(params.get_use_mp() ) {
		hmcout << "\t" << iter0 << "\t" << iter1;
	}
	if(meta::get_use_rectangles(params) ) {
		//print rectangle value
		hmcout << "\t" << observables.rectangles;
	}
	hmcout << std::endl;
	hmcout.close();

	//print to screen
	print_hmcobservables();
}

void hmcExecutable::print_hmcobservables()
{
	using namespace std;
	//short version of output, all obs are collected in the output file anyways...
	logger.info() << "\tHMC [OBS]:\t" << iteration << setw(8) << setfill(' ') << "\t" << setprecision(15) << observables.plaq << "\t" << observables.poly.re << "\t" << observables.poly.im;
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
