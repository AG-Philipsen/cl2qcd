#include "heatbathExecutable.h"

inline heatbathExecutable::heatbathExecutable(int argc, const char* argv[]) :
		generationExecutable(argc, argv)
{
	initializationTimer.reset();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	meta::print_info_heatbath(ownName, parameters);
	writeHeatbathLogfile();
	setIterationParameters();
	initializationTimer.add();
}

inline void heatbathExecutable::writeHeatbathLogfile()
{
	outputToFile.open(filenameForHeatbathLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_heatbath(ownName, &outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForHeatbathLogfile);
	}
}

void heatbathExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps = parameters.get_heatbathsteps();
	overrelaxSteps = parameters.get_overrelaxsteps();
}

void heatbathExecutable::performThermalization()
{
	logger.info() << "Start thermalization";
	int iteration = 0;
	writeGaugeObservablesToScreen(iteration);
	for (; iteration < thermalizationSteps; iteration++)
	{
		physics::algorithms::heatbath(*gaugefield, *prng);
	}
	logger.info() << "thermalization done";
}

inline void heatbathExecutable::performHeatbathAndMeasurements()
{
	logger.info() << "Start heatbath";
	for (int iteration = 0; iteration < generationSteps; iteration++)
	{
		physics::algorithms::heatbath(*gaugefield, *prng, overrelaxSteps);
		measureGaugeObservables(iteration);
		saveGaugefield(iteration);
	}
	logger.info() << "heatbath done";
}

int main(int argc, const char* argv[])
{
	try {
		heatbathExecutable heatbathInstance(argc, argv);
		heatbathInstance.generateConfigurations();
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
