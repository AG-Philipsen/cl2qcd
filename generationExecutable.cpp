#include "generationExecutable.h"

generationExecutable::generationExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
{
	initializationTimer.reset();
	gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	initializationTimer.add();
}

void generationExecutable::setIterationParameters()
{
	//NOTE: this is 0 in case of cold or hot start
	iteration 				= gaugefield->get_parameters_source().trajectorynr_source;
	thermalizationSteps 	= iteration + parameters.get_thermalizationsteps();
	generationSteps 		= thermalizationSteps;
	writeFrequency 			= parameters.get_writefrequency();
	saveFrequency 			= parameters.get_savefrequency();
}

void generationExecutable::measureGaugeObservables()
{
  gaugeObservablesInstance.measurePlaqAndPoly(*gaugefield, iteration, parameters);
  if ( parameters.get_measure_transportcoefficient_kappa() ) {
    measureTransportcoefficientKappa();
  }
}

void generationExecutable::saveGaugefield()
{
        gaugefield->save(iteration+1);
	if (((saveFrequency != 0) && ((iteration + 1) % saveFrequency) == 0)) {
		gaugefield->saveToSpecificFile(iteration + 1);
	}
}

void generationExecutable::savePrng()
{
	prng->save(iteration + 1);
	if (((saveFrequency != 0) && ((iteration + 1) % saveFrequency) == 0)) {
		prng->saveToSpecificFile(iteration + 1);
	}
}

void generationExecutable::measureTransportcoefficientKappa()
{
	double kappa = 0;
	kappa = physics::algorithms::kappa_clover(*gaugefield, parameters.get_beta());
	writeTransportcoefficientKappaToFile(kappa, "kappa_clover.dat");
}

void generationExecutable::writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa)
{
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << iteration << "\t" << kappa << std::endl;
}

void generationExecutable::writeTransportcoefficientKappaToFile(hmc_float kappa, std::string filename)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if ( outputToFile.is_open() ) {
		writeTransportcoefficientKappaToFileUsingOpenOutputStream(kappa);
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
	measureGaugeObservables();
	for (; iteration < thermalizationSteps; iteration++)
	 {
		thermalizeAccordingToSpecificAlgorithm();
	}
	logger.info() << "...thermalization done";
}

void generationExecutable::generate()
{
	logger.info() << "Start generation of configurations...";
	for (; iteration < generationSteps; iteration++)
	 {
		generateAccordingToSpecificAlgorithm();
		performOnlineMeasurements();
		saveGaugefield();
		savePrng();
	}
	logger.info() << "...generation done";
}

void generationExecutable::performOnlineMeasurements()
{
	if( ( (iteration + 1) % writeFrequency ) == 0 ) {
		measureGaugeObservables();
	}
}

