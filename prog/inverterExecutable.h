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
	 * Writes Logfile.
	 */
	inverterExecutable(int argc, const char* argv[]);

protected:
	const std::string 	filenameForInverterLogfile = "inverter.log";
	std::string 		filenameForTwoFlavourDoubletChiralCondensateData;
	std::string 		filenameForTwoFlavourDoubletCorrelatorData;

	void writeInverterLogfile();

	void printParametersToScreenAndFile();

	void measureTwoFlavourDoubletCorrelatorsOnGaugefield();

	void measureTwoFlavourDoubletChiralCondensateOnGaugefield();

	/**
	 * Performs measurements of fermionic observables on possibly multiple gaugefield configurations.
	 */
	void performApplicationSpecificMeasurements();
};

#endif /* _INVERTERH_ */

