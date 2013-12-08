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
	gaugeObservablesInstance.measureGaugeObservables(*gaugefield, iteration, parameters);
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
	  gaugeObservablesInstance.measureGaugeObservables(*gaugefield, iteration, parameters);
	}
}

