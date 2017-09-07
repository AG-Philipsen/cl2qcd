/**
 * Copyright 2015 Christopher Pinke
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

# pragma once

#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"

namespace hardware
{
	class HardwareParametersInterface
	{
	public:
		HardwareParametersInterface() {};
		virtual ~HardwareParametersInterface() {};
		virtual int getNs() const = 0;
		virtual int getNt() const = 0;
		virtual int getSpatialLatticeVolume() const = 0;
		virtual int getLatticeVolume() const = 0;
		virtual bool useGpu() const = 0;
		virtual bool useCpu() const = 0;
		virtual int getMaximalNumberOfDevices() const = 0;
		virtual std::vector<int> getSelectedDevices() const = 0;
		virtual bool splitCpu() const = 0;
		virtual bool enableProfiling() const = 0;
		virtual bool disableOpenCLCompilerOptimizations() const = 0;
		virtual bool useSameRandomNumbers() const = 0;
		virtual bool useEvenOddPreconditioning() const = 0;
	};
}
