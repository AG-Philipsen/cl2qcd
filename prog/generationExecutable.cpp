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
