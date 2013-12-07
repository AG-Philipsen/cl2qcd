/** @file
 *
 * Everything required by heatbath's main()
 */

#ifndef HEATBATHEXECUTABLE_H_
#define HEATBATHEXECUTABLE_H_

#include "generationExecutable.h"

#include "meta/util.hpp"
#include "physics/prng.hpp"
#include "physics/lattices/gaugefield.hpp"
#include "physics/algorithms/heatbath.hpp"

class heatbathExecutable: public generationExecutable
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
	int heatbathSteps;
	int overrelaxSteps;

	void writeHeatbathLogfile();

	void setIterationParameters();

	void performThermalization();

	void performHeatbathAndMeasurements();
};

#endif /* HEATBATHEXECUTABLE_H_ */
