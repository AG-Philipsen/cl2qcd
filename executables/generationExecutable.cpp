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
	iteration           = gaugefield->get_parameters_source().trajectorynr_source;
	thermalizationSteps = iteration + parameters.get_thermalizationsteps();
	generationSteps     = thermalizationSteps; //this is temp.: in each child it is incremented by the nr of tr.
	writeFrequency      = parameters.get_writefrequency();
	saveFrequency       = parameters.get_savefrequency();
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
	physics::observables::measureGaugeObservablesAndWriteToFile(gaugefield, iteration);
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
          physics::observables::measureGaugeObservablesAndWriteToFile(gaugefield, iteration);
	}
}

