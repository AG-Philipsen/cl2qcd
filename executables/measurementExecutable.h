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


/*
 * @file
 * Declaration of the measurementExecutable class.
 * This class provides features for measurements on (multiple) gauge configurations.
 */

#ifndef MEASUREMENTEXECUTABLE_H_
#define MEASUREMENTEXECUTABLE_H_

#include "generalExecutable.h"

class measurementExecutable : public generalExecutable
{
public:
	measurementExecutable(int argc, const char* argv[]);

	void performMeasurements();

protected:
	std::string currentConfigurationName;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;
	int iteration;

	void setIterationVariables();

	void checkStartconditions();

	void initializeGaugefieldAccordingToIterationVariable();

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();

	void initializeGaugefield();

	void performMeasurementsForSpecificIteration();

	virtual void performApplicationSpecificMeasurements() {};
};


#endif /* MEASUREMENTEXECUTABLE_H_ */
