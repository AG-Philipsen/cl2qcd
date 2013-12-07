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

	void setIterationVariables();

	void initializeGaugefieldAccordingToIterationVariable(int interation);

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();

	void initializeGaugefield(int interation);

	void performMeasurementsForSpecificIteration(int interation);

	virtual void performApplicationSpecificMeasurements() {};
};


#endif /* MEASUREMENTEXECUTABLE_H_ */
