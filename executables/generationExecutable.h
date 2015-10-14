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


/**
 * @file
 * Declaration of the generationExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to certain algorithms.
 **/

#include "./generalExecutable.h"

class generationExecutable : public generalExecutable
{
public:
	void generateConfigurations();
	virtual ~generationExecutable(){};

protected:
	//Protected since it makes no sense to allow the user to instatiate this class
	generationExecutable(int argc, const char* argv[], std::string parameterSet = "all parameters");
	
	int writeFrequency;
	int saveFrequency;
    int savePointFrequency;
	int thermalizationSteps;
	int generationSteps;
	int iteration;
// 	std::string filenameForGaugeobservables;

	/**
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 **/
	void setIterationParameters();

	/**
	 * Saves current gaugefield configuration to disk if current iteration
	 * is a multiple of the inputparameter save_frequency or if it is the
	 * last iteration.
	 **/
	void saveGaugefield();

	/**
	 * Saves current prng configuration to disk if current iteration
	 * is a multiple of the inputparameter save_frequency or if it is the
	 * last iteration.
	 **/
	void savePrng();

	/**
	 * Performs thermalization of the physical system according to the algorithm
	 * specified in the "thermalizeAccordingToSpecificAlgorithm" function.
	 **/
	void thermalize();

	/**
	 * Generates a new gauge configuration according to the algorithm
	 * specified in the "generateAccordingToSpecificAlgorithm" function.
	 * Performs measurements on this configuration according to the
	 * function "performOnlineMeasurements".
	 **/
	void generate();

	void virtual thermalizeAccordingToSpecificAlgorithm() = 0;

	void virtual generateAccordingToSpecificAlgorithm() = 0;

	/**
	 * Measurements to be performed after each step of configuration generation.
	 * By default this measures the gauge observables.
	 * This function can be replaced by a more specific one for each specific algorithm.
	 */
	void virtual performOnlineMeasurements();
};


