/*
 * @file
 * Declaration of the measurementExecutable class.
 * This class provides features for measurements on (multiple) gauge configurations.
 */

#ifndef MEASUREMENTEXECUTABLE_H_
#define MEASUREMENTEXECUTABLE_H_

#include "./executables/generalExecutable.h"

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
