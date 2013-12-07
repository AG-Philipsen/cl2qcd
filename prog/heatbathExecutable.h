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
	heatbathExecutable(int argc, const char* argv[]);

private:
	const std::string filenameForHeatbathLogfile = "heatbath.log";
	int heatbathSteps;
	int overrelaxSteps;

	void thermalizeAccordingToSpecificAlgorithm();

	void generateAccordingToSpecificAlgorithm();

	void writeHeatbathLogfile();

	void setIterationParameters();
};

#endif /* HEATBATHEXECUTABLE_H_ */
