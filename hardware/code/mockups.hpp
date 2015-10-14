/*
 * Copyright 2015 Christopher Pinke, Francesca Cuteri
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

#include "../hardwareParameters.hpp"

class HardwareParametersMockup : public hardware::HardwareParametersInterface
{
public:
	HardwareParametersMockup(const int nsIn, const int ntIn) : ns(nsIn), nt(ntIn) {};
	~HardwareParametersMockup() {};
	virtual int getNs() const override
	{
		return ns;
	}
	virtual int getNt() const override
	{
		return nt;
	}
	virtual bool disableOpenCLCompilerOptimizations() const override
	{
		return false;
	}
	virtual bool useGpu() const override
	{
		return false;
	}
	virtual bool useCpu() const override
	{
		return true;
	}
	virtual int getMaximalNumberOfDevices() const override
	{
		return 1;
	}
	virtual std::vector<int> getSelectedDevices() const override
	{
		return std::vector<int>{0};
	}
	virtual bool splitCpu() const override
	{
		return false;
	}
	virtual bool enableProfiling() const override
	{
		return false;
	}
	virtual bool useSameRandomNumbers() const override
	{
		return false;
	}
	virtual bool useEvenOddPreconditioning() const override
	{
		return false;
	}
	virtual int getSpatialLatticeVolume() const override
	{
		return getNs() * getNs() * getNs();
	}
	virtual int getLatticeVolume() const override
	{
		return getNs() * getNs() * getNs() * getNt();
	}
private:
	int ns, nt;
};

