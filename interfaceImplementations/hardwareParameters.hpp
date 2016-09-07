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

#include "../hardware/hardwareParameters.hpp"
#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"

namespace hardware
{
	class HardwareParametersImplementation final : public HardwareParametersInterface
	{
	public:
		HardwareParametersImplementation( const meta::Inputparameters * parametersIn) : fullParameters(parametersIn) {}
		~HardwareParametersImplementation() {};
		virtual int getNs() const override
		{
			return fullParameters->get_nspace();
		}
		virtual int getNt() const override
		{
			return fullParameters->get_ntime();
		}
		virtual bool disableOpenCLCompilerOptimizations() const override
		{
			return fullParameters->is_ocl_compiler_opt_disabled();
		}
		virtual bool useGpu() const override
		{
			return fullParameters->get_use_gpu();
		}
		virtual bool useCpu() const override
		{
			return fullParameters->get_use_cpu();
		}
		virtual int getMaximalNumberOfDevices() const override
		{
			return fullParameters->get_device_count();
		}
		virtual std::vector<int> getSelectedDevices() const override
		{
			return fullParameters->get_selected_devices();
		}
		virtual bool splitCpu() const override
		{
			return fullParameters->get_split_cpu();
		}
		virtual bool enableProfiling() const override
		{
			return fullParameters->get_enable_profiling();
		}
		virtual bool useSameRandomNumbers() const override
		{
			return fullParameters->get_use_same_rnd_numbers();
		}
		virtual bool useEvenOddPreconditioning() const override
		{
			return fullParameters->get_use_eo();
		}
		virtual int getSpatialLatticeVolume() const override
		{
			return meta::get_volspace( *fullParameters );
		}
		virtual int getLatticeVolume() const override
		{
			return meta::get_vol4d( * fullParameters );
		}
	private:
		const meta::Inputparameters * fullParameters;
	};
}
