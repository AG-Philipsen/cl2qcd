#include "generationExecutable.h"

generationExecutable::generationExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
{
	initializationTimer.reset();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	initializationTimer.add();
}

void generationExecutable::setIterationParameters()
{
	thermalizationSteps = parameters.get_thermalizationsteps();
	writeFrequency = parameters.get_writefrequency();
	saveFrequency = parameters.get_savefrequency();
}

void generationExecutable::writeGaugeObservablesToFile(int& iteration)
{
	if (((iteration + 1) % writeFrequency) == 0) {
		filenameForGaugeobservables = meta::get_gauge_obs_file_name(parameters,
				"");
		print_gaugeobservables(*gaugefield, iteration,
				filenameForGaugeobservables);
	}
}

void generationExecutable::writeGaugeObservablesToScreen(int& iteration)
{
	if (parameters.get_print_to_screen() || (iteration == 0)) {
		print_gaugeobservables(*gaugefield, iteration);
	}
}

void generationExecutable::writeGaugeObservablesToScreenAndFile(int iteration) {
	writeGaugeObservablesToScreen(iteration);
	writeGaugeObservablesToFile(iteration);
}

void generationExecutable::measureGaugeObservables(int& iteration)
{
	writeGaugeObservablesToScreenAndFile(iteration);
	if ( parameters.get_measure_transportcoefficient_kappa() ) {
		measureTransportcoefficientKappa(iteration);
	}
}

void generationExecutable::saveGaugefield(int iteration)
{
	if (((saveFrequency != 0) && ((iteration + 1) % saveFrequency) == 0)) {
		gaugefield->save(iteration + 1);
	}
	if (iteration == generationSteps - 1) {
		gaugefield->save("conf.save", iteration + 1);
	}
}

void generationExecutable::measureTransportcoefficientKappa(int iteration)
{
	double kappa = 0;
	kappa = physics::algorithms::kappa_clover(*gaugefield, parameters.get_beta());
	writeTransportcoefficientKappaToFile(kappa, iteration, "kappa_clover.dat");
}

void generationExecutable::writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa, int iteration)
{
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << iteration << "\t" << kappa << std::endl;
}

void generationExecutable::writeTransportcoefficientKappaToFile(hmc_float kappa, int iteration, std::string filename)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if ( outputToFile.is_open() ) {
		writeTransportcoefficientKappaToFileUsingOpenOutputStream(kappa, iteration);
		outputToFile.close();
	} else {
		logger.warn() << "Could not open " << filename;
		File_Exception(filename.c_str());
	}
}

void generationExecutable::generateConfigurations()
{
	performanceTimer.reset();
	thermalize();
	generate();
	performanceTimer.add();
}

void generationExecutable::thermalize()
{
	logger.info() << "Start thermalization...";
	int iteration = 0;
	writeGaugeObservablesToScreen(iteration);
	for (; iteration < thermalizationSteps; iteration++)
	 {
		thermalizeAccordingToSpecificAlgorithm();
	}
	logger.info() << "...thermalization done";
}

void generationExecutable::generate()
{
	logger.info() << "Start generation of configurations...";
	for (int iteration = 0; iteration < generationSteps; iteration++)
	 {
		generateAccordingToSpecificAlgorithm();
		performOnlineMeasurements(iteration);
		saveGaugefield(iteration);
	}
	logger.info() << "...generation done";
}


