/** @file
 *
 * Everything required by heatbath's main()
 */

#ifndef HEATBATHEXECUTABLE_H_
#define HEATBATHEXECUTABLE_H_

#include "generalExecutable.h"

#include "meta/util.hpp"
#include "physics/prng.hpp"
#include "physics/lattices/gaugefield.hpp"
#include "physics/algorithms/heatbath.hpp"
#include "physics/algorithms/kappa_clover.hpp"

class heatbathExecutable: public generalExecutable
{
public:
	/**
	 * Constructor.
	 * Initializes generalExecutable object, Random Number Generator and Iteration variables.
	 */
	heatbathExecutable(int argc, const char* argv[]);

	/**
	 * Performs heatbath algorithms (including thermalization) on gaugefield object and
	 * measures of gauge observables on each iteration (not during thermalization).
	 */
	void performHeatbathAndMeasureGaugeObservables();

private:
	const std::string 	filenameForHeatbathLogfile 		= "heatbath.log";
	std::string filenameForGaugeobservables;
	int thermalizationSteps;
	int heatbathSteps;
	int overrelaxSteps;
	int writeFrequency;
	int saveFrequency;

	void writeHeatbathLogfile();

	void setIterationParameters();

	void performThermalization();

	void writeGaugeObservablesToFile(int& iteration);

	void writeGaugeObservablesToScreen(int& iteration);

	void writeGaugeObservablesToScreenAndFile(int iteration);

	void saveGaugefield(int iteration);

	void measureGaugeObservables(int& iteration);

	void performHeatbathAndMeasurements();

	void measureTransportcoefficientKappa(int iteration);

	void writeTransportcoefficientKappaToFile(hmc_float kappa, int iteration, std::string filename);

	void writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa, int iteration);
};

#endif /* HEATBATHEXECUTABLE_H_ */
