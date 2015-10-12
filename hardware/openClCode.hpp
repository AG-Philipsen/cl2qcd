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

#include "code/gaugefield.hpp"
#include "code/prng.hpp"
#include "code/real.hpp"
#include "code/complex.hpp"
#include "code/spinors.hpp"
#include "code/spinors_staggered.hpp"
#include "code/fermions.hpp"
#include "code/fermions_staggered.hpp"
#include "code/correlator.hpp"
#include "code/correlator_staggered.hpp"
#include "code/heatbath.hpp"
#include "code/kappa.hpp"
#include "code/gaugemomentum.hpp"
#include "code/molecular_dynamics.hpp"
#include "code/buffer.hpp"

#include <memory>

//Todo: make device arg constant

namespace hardware
{
	class OpenClCode
	{
	public:
		virtual ~OpenClCode(){};
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Gaugefield> getCode_gaugefield(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Complex> getCode_complex(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Prng> getCode_PRNG(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Spinors> getCode_Spinors(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Fermions> getCode_Fermions(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Gaugemomentum> getCode_Gaugemomentum(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Molecular_Dynamics> getCode_Molecular_Dynamics(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Correlator> getCode_Correlator(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Heatbath> getCode_Heatbath(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Kappa> getCode_Kappa(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Buffer> getCode_Buffer(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Spinors_staggered> getCode_Spinors_staggered(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Correlator_staggered> getCode_Correlator_staggered(hardware::Device *) const = 0;
		virtual std::unique_ptr<const hardware::code::Fermions_staggered> getCode_Fermions_staggered(hardware::Device *) const = 0;
	};
}

#include "../meta/inputparameters.hpp"

namespace hardware
{
	class OpenClCode_fromMetaInputparameters final : public OpenClCode
	{
	public:
		OpenClCode_fromMetaInputparameters( const meta::Inputparameters & parametersIn ) : parameters(parametersIn), kernelParameters(nullptr)
		{
			kernelParameters = new hardware::code::OpenClKernelParametersImplementation( parametersIn );
		};
		~OpenClCode_fromMetaInputparameters()
		{
			if (kernelParameters)
			{
				delete kernelParameters;
			}
		}
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Real>( new hardware::code::Real{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugefield> getCode_gaugefield(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugefield>( new hardware::code::Gaugefield{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Prng> getCode_PRNG(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Prng>( new hardware::code::Prng{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Complex> getCode_complex(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Complex>( new hardware::code::Complex{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors> getCode_Spinors(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors>( new hardware::code::Spinors{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions> getCode_Fermions(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions>( new hardware::code::Fermions{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugemomentum> getCode_Gaugemomentum(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugemomentum>( new hardware::code::Gaugemomentum{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Molecular_Dynamics> getCode_Molecular_Dynamics(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Molecular_Dynamics>( new hardware::code::Molecular_Dynamics{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator> getCode_Correlator(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator>( new hardware::code::Correlator{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Heatbath> getCode_Heatbath(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Heatbath>( new hardware::code::Heatbath{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Kappa> getCode_Kappa(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Kappa>( new hardware::code::Kappa{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Buffer> getCode_Buffer(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Buffer>( new hardware::code::Buffer{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors_staggered> getCode_Spinors_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors_staggered>( new hardware::code::Spinors_staggered{parameters, *kernelParameters,deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator_staggered> getCode_Correlator_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator_staggered>( new hardware::code::Correlator_staggered{parameters, *kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions_staggered> getCode_Fermions_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions_staggered>( new hardware::code::Fermions_staggered{parameters, *kernelParameters, deviceIn} ) ;
		}
	private:
		const hardware::code::OpenClKernelParametersInterface * kernelParameters;
		const meta::Inputparameters & parameters;
	};
}
