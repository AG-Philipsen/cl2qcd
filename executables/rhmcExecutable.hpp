/*
 * Copyright (c) 2014,2018 Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @file
 * Declaration of the rhmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Rational Hybrid Monte Carlo (RHMC) algorithm.
 */

#ifndef RHMCEXECUTABLE_H_
#define RHMCEXECUTABLE_H_

#include "../physics/algorithms/rhmc.hpp"
#include "../physics/algorithms/rational_approximation.hpp"
#include "../physics/observables/staggeredChiralCondensate.hpp"
#include <cmath>
#include "generationExecutable.hpp"

class rhmcExecutable : public generationExecutable {
public:
	rhmcExecutable(int argc, const char* argv[]);
	~rhmcExecutable();

protected:
	double acceptanceRate = 0;
	hmc_observables observables;

	/*
	 * Sets member variables that control the iterations during
	 * the generation of gaugefield configurations.
	 */
	void setIterationParameters();

	void printParametersToScreenAndFile();

	void writeRhmcLogfile();

	void thermalizeAccordingToSpecificAlgorithm() override;

	void generateAccordingToSpecificAlgorithm() override;

	/**
	 * Measures RHMC related observables
	 */
	void performOnlineMeasurements() override;

	void printRhmcObservables(const std::string& filename);

	void printRhmcObservablesToFile(const std::string& filename);

	void printRhmcObservablesToScreen();

private:
	physics::algorithms::Rational_Approximation *approx_hb, *approx_md, *approx_met;
	void checkRhmcParameters(const meta::Inputparameters& p);

};

#endif /* RHMCEXECUTABLE_H_ */
