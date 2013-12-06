/** @file
 *
 * Everything required by inverter's main()
 */
#ifndef _INVERTERH_
#define _INVERTERH_

#include "generalExecutable.h"
#include "physics/lattices/spinorfield.hpp"
#include "physics/sources.hpp"
#include "physics/algorithms/flavour_doublet.hpp"
#include "physics/algorithms/inversion.hpp"

/**
 * Inverter executable which measures fermionic observables on gaugefield configurations.
 */
class inverterExecutable : public multipleConfigurationExecutable
{
public:
	/**
	 * Constructor.
	 * Initializes generalExecutable object, Random Number Generator and Iteration variables.
	 */
	inverterExecutable(int argc, const char* argv[]);

	/**
	 * Performs measurements of fermionic observables on possibly multiple gaugefield configurations.
	 */
	void performMeasurements();

protected:
	const std::string 	filenameForInverterLogfile 		= "inverter.log";
	std::string 		filenameForTwoFlavourDoubletChiralCondensateData;
	std::string 		filenameForTwoFlavourDoubletCorrelatorData;
	std::string 		currentConfigurationName;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;

	void setIterationVariables();

	void initializeGaugefieldAccordingToIterationVariable(int interation);

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();

	void initializeGaugefield(int interation);

	void performMeasurementsForSpecificIteration(int interation);

	void writeInverterLogfile();

	void measureTwoFlavourDoubletCorrelatorsOnGaugefield();

	void measureTwoFlavourDoubletChiralCondensateOnGaugefield();

	void measureFermionicObservablesOnGaugefield();
};

#endif /* _INVERTERH_ */

