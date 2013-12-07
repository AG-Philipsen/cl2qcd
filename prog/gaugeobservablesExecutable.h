/** @file
 *
 * Everything required by gaugeobservable's main()
 */

#ifndef GAUGEOBSERVABLESEXECUTABLE_H_
#define GAUGEOBSERVABLESEXECUTABLE_H_

#include "generalExecutable.h"

class gaugeobservablesExecutable : public measurementExecutable
{
public:
	gaugeobservablesExecutable(int argc, const char* argv[]);

protected:
	std::string filenameForGaugeobservables;
	const std::string 	filenameForGaugeobservablesLogfile = "gaugeobservables.log";

	void writeGaugeobservablesLogfile();

	void printParametersToScreenAndFile();

	/**
	 * Performs measurements of gauge observables on possibly multiple gaugefield configurations.
	 */
	void performApplicationSpecificMeasurements();
};

#endif /* GAUGEOBSERVABLESEXECUTABLE_H_ */
