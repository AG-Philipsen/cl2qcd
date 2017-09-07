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

generationExecutable::generationExecutable(int argc, const char* argv[], std::string parameterSet) : generalExecutable(argc, argv, parameterSet)
{
	initializationTimer.reset();
	gaugefield = new physics::lattices::Gaugefield(*system, &(interfacesHandler->getInterface<physics::lattices::Gaugefield>()), *prng);
	initializationTimer.add();
}

void generationExecutable::setIterationParameters()
{
	//NOTE: this is 0 in case of cold or hot start
	iteration           = gaugefield->get_trajectoryNumberAtInit();
	thermalizationSteps = iteration + parameters.get_thermalizationsteps();
	generationSteps     = thermalizationSteps; //this is temp.: in each child it is incremented by the nr of tr.
	writeFrequency      = parameters.get_writefrequency();
	saveFrequency       = parameters.get_savefrequency();
    savePointFrequency  = parameters.get_savepointfrequency();
}

void generationExecutable::saveGaugefield()
{
    if (((savePointFrequency != 0) && ((iteration + 1) % savePointFrequency) == 0)) {
        gaugefield->save(iteration+1); //Here the number is that written in the lime file as metadata, and it is iteration+1 to be able to continue later at the right tr.
    }
    if (((saveFrequency != 0) && ((iteration + 1) % saveFrequency) == 0)) {
		gaugefield->saveToSpecificFile(iteration + 1);
	}
}

void generationExecutable::savePrng()
{
    if (((savePointFrequency != 0) && ((iteration + 1) % savePointFrequency) == 0)) {
        prng->save();
    }
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
	logger.info() << "Start thermalization (" << parameters.get_thermalizationsteps() << " tr.)...";
	physics::observables::measureGaugeObservablesAndWriteToFile(gaugefield, iteration, interfacesHandler->getGaugeObservablesParametersInterface());
	//With this try and catch the warning is printed only if the user wants to make thermalization steps, not always, but this makes sense
	try
	{
	      for (; iteration < thermalizationSteps; iteration++)
	      {
		    thermalizeAccordingToSpecificAlgorithm();
	      }
	}
	catch(Print_Error_Message& exception){
   	      logger.warn() << "The thermalization is not yet implemented!  It is just skipped.";
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
          physics::observables::measureGaugeObservablesAndWriteToFile(gaugefield, iteration, interfacesHandler->getGaugeObservablesParametersInterface());
	}
}

