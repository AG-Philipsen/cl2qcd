/*
 * @file
 * Declaration of the generationExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to certain algorithms.
 */

#include "generalExecutable.h"
#include "physics/algorithms/kappa_clover.hpp"

class generationExecutable : public generalExecutable
{
public:
	generationExecutable(int argc, const char* argv[]);

	void generateConfigurations();

protected:
	int writeFrequency;
	int saveFrequency;
	int thermalizationSteps;
	int generationSteps;
	std::string filenameForGaugeobservables;

	void setIterationParameters();

	void writeGaugeObservablesToFile(int& iteration);

	void writeGaugeObservablesToScreen(int& iteration);

	void writeGaugeObservablesToScreenAndFile(int iteration);

	void measureGaugeObservables(int& iteration);

	void saveGaugefield(int iteration);

	void measureTransportcoefficientKappa(int iteration);

	void writeTransportcoefficientKappaToFile(hmc_float kappa, int iteration, std::string filename);

	void writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa, int iteration);

	void virtual thermalize() {};

	void virtual generate() {};
};


