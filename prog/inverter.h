/** @file
 *
 * Everything required by inverter's main()
 */
#ifndef _INVERTERH_
#define _INVERTERH_

#include "generalExecutable.h"
#include "physics/lattices/gaugefield.hpp"
#include "physics/lattices/spinorfield.hpp"
#include "physics/sources.hpp"
#include "physics/algorithms/flavour_doublet.hpp"
#include "physics/algorithms/inversion.hpp"

class inverterExecutable : public generalExecutable
{
public:
	inverterExecutable(int argc, const char* argv[]);

	~inverterExecutable();

	void performMeasurements();

protected:
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	const std::string 	filenameForCurrentPrngState 	= "prng.inverter.save";
	const std::string 	filenameForInverterLogfile 		= "inverter.log";
	const std::string 	filenameForProfilingData 		= std::string(ownName) + std::string("_profiling_data");
	std::string 		filenameForTwoFlavourDoubletChiralCondensateData;
	std::string 		filenameForTwoFlavourDoubletCorrelatorData;
	std::string 		currentConfigurationName;
	std::fstream 		outputStreamForProfilingData;
	int iterationStart;
	int iterationEnd;
	int iterationIncrement;
	int iteration;
	usetimer solverTimer;

	void setIterationVariables();

	void saveCurrentPrngStateToFile();

	void initializeGaugefieldAccordingToIterationVariable();

	void initializeGaugefieldAccordingToConfigurationGivenInSourcefileParameter();

	void initializeGaugefield();

	void performMeasurementsForSpecificIteration();

	void writeInverterLogfile();

	void writeProfilingDataToScreenAndFile();

	void getSolverStatistics(int& totalSolverCalls, uint64_t& totalSolverTime, uint64_t& averageSolverTime);

	void writeProfilingDataToScreen(uint64_t totalSolverTime, int totalSolverCalls, uint64_t averageSolverTime);

	void writeProfilingDataToFile(uint64_t totalSolverTime,	int totalSolverCalls, uint64_t averageSolverTime);

	void writeProfilingDataToFileUsingOpenOutputStream(uint64_t totalSolverTime, int totalSolverCalls, uint64_t averageSolverTime);

	void measureTwoFlavourDoubletCorrelatorsOnGaugefield();

	void measureTwoFlavourDoubletChiralCondensateOnGaugefield();

	void measureFermionicObservablesOnGaugefield();
};

#endif /* _INVERTERH_ */

