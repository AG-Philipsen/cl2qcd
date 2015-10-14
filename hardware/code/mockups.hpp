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

namespace hardware
{
class HardwareParametersMockup : public HardwareParametersInterface
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
}
#include "../openClCode.hpp"

namespace hardware
{
	class OpenClCodeMockup : public OpenClCode
	{
	public:
		OpenClCodeMockup() : kernelParameters(nullptr)
		{};
		~OpenClCodeMockup()
		{
			if (kernelParameters)
			{
				delete kernelParameters;
			}
		}
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Real>( new hardware::code::Real{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugefield> getCode_gaugefield(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugefield>( new hardware::code::Gaugefield{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Prng> getCode_PRNG(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Prng>( new hardware::code::Prng{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Complex> getCode_complex(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Complex>( new hardware::code::Complex{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors> getCode_Spinors(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors>( new hardware::code::Spinors{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions> getCode_Fermions(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions>( new hardware::code::Fermions{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugemomentum> getCode_Gaugemomentum(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugemomentum>( new hardware::code::Gaugemomentum{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Molecular_Dynamics> getCode_Molecular_Dynamics(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Molecular_Dynamics>( new hardware::code::Molecular_Dynamics{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator> getCode_Correlator(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator>( new hardware::code::Correlator{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Heatbath> getCode_Heatbath(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Heatbath>( new hardware::code::Heatbath{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Kappa> getCode_Kappa(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Kappa>( new hardware::code::Kappa{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Buffer> getCode_Buffer(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Buffer>( new hardware::code::Buffer{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors_staggered> getCode_Spinors_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors_staggered>( new hardware::code::Spinors_staggered{*kernelParameters,deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator_staggered> getCode_Correlator_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator_staggered>( new hardware::code::Correlator_staggered{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions_staggered> getCode_Fermions_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions_staggered>( new hardware::code::Fermions_staggered{*kernelParameters, deviceIn} ) ;
		}
	private:
		const hardware::code::OpenClKernelParametersInterface * kernelParameters;
	};
}
