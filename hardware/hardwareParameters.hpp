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

class HardwareParametersInterface
{
protected:
	HardwareParametersInterface() {};
	virtual ~HardwareParametersInterface() {};
	virtual bool useGpu() const = 0;
	virtual bool useCpu() const = 0;
	virtual int getMaximalNumberOfDevices() const = 0;
	virtual std::vector<int> getSelectedDevices() const = 0;
	virtual bool splitCpu() const = 0;
	virtual bool enableProfiling() const = 0;
};

class HardwareParameters : public HardwareParametersInterface
{
public:
	HardwareParameters( meta::Inputparameters * parametersIn) : fullParameters(parametersIn) {}
	~HardwareParameters() {};
	bool useGpu() const
	{
		return fullParameters->get_use_gpu();
	}
	virtual bool useCpu() const
	{
		return fullParameters->get_use_cpu();
	}
	virtual int getMaximalNumberOfDevices() const
	{
		return fullParameters->get_device_count();
	}
	virtual std::vector<int> getSelectedDevices() const
	{
		return fullParameters->get_selected_devices();
	}
	virtual bool splitCpu() const
	{
		return fullParameters->get_split_cpu();
	}
	virtual bool enableProfiling() const
	{
		return fullParameters->get_enable_profiling();
	}
private:
	meta::Inputparameters * fullParameters;
};
