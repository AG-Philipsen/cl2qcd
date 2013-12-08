/** @file
 *
 * Everything required by heatbath's main()
 */

#ifndef HEATBATHEXECUTABLE_H_
#define HEATBATHEXECUTABLE_H_

#include "generationExecutable.h"
#include "physics/algorithms/heatbath.hpp"

class heatbathExecutable: public generationExecutable
{
public:
	heatbathExecutable(int argc, const char* argv[]);

private:
	const std::string filenameForHeatbathLogfile = "heatbath.log";
	int overrelaxSteps;

	/*
	 * Thermalize the system using the heatbath algorithm.
	 */
	void thermalizeAccordingToSpecificAlgorithm();

	/*
	 * Generate configurations using the heatbath algorithm.
	 */
	void generateAccordingToSpecificAlgorithm();

	void writeHeatbathLogfile();

	void setIterationParameters();

	void performOnlineMeasurements(int iteration);
};

#endif /* HEATBATHEXECUTABLE_H_ */
