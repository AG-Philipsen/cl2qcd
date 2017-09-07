/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * @file
 * Declaration of the hmcExecutable class.
 * This class provides features for the generation of gauge configurations
 * according to the Hybrid Monte Carlo (HMC) algorithm.
 */

#ifndef HMCEXECUTABLE_H_
#define HMCEXECUTABLE_H_

#include "generationExecutable.h"
#include "../physics/algorithms/hmc.hpp"
#include <cmath>

class hmcExecutable : public generationExecutable {
    public:
        hmcExecutable(int argc, const char* argv[]);
        ~hmcExecutable();
    protected:
        // 	const std::string filenameForHmcLogfile = "hmc.log";
        double acceptanceRate = 0;
        hmc_observables observables;

        /*
         * Sets member variables that control the iterations during
         * the generation of gaugefield configurations.
         */
        void setIterationParameters();

        void printParametersToScreenAndFile();

        void writeHmcLogfile();

        void thermalizeAccordingToSpecificAlgorithm() override;

        void generateAccordingToSpecificAlgorithm() override;

        /**
         * Measures HMC related observables
         */
        void performOnlineMeasurements() override;

        void printHmcObservables(const std::string& filename);

        void printHmcObservablesToFile(const std::string& filename);

        void printHmcObservablesToScreen();
};

#endif /* HMCEXECUTABLE_H_ */
