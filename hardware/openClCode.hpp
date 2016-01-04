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

namespace hardware
{
	class OpenClCode
	{
	public:
		OpenClCode( const hardware::code::OpenClKernelParametersInterface & kernelParametersIn ) : kernelParameters(&kernelParametersIn)
		{};
		~OpenClCode()
		{}
		std::unique_ptr<hardware::code::Real> getCode_real(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Real>( new hardware::code::Real{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Gaugefield> getCode_gaugefield(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Gaugefield>( new hardware::code::Gaugefield{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Prng> getCode_PRNG(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Prng>( new hardware::code::Prng{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Complex> getCode_complex(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Complex>( new hardware::code::Complex{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Spinors> getCode_Spinors(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Spinors>( new hardware::code::Spinors{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Fermions> getCode_Fermions(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Fermions>( new hardware::code::Fermions{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Gaugemomentum> getCode_Gaugemomentum(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Gaugemomentum>( new hardware::code::Gaugemomentum{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Molecular_Dynamics> getCode_Molecular_Dynamics(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Molecular_Dynamics>( new hardware::code::Molecular_Dynamics{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Correlator> getCode_Correlator(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Correlator>( new hardware::code::Correlator{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Heatbath> getCode_Heatbath(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Heatbath>( new hardware::code::Heatbath{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Kappa> getCode_Kappa(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Kappa>( new hardware::code::Kappa{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Buffer> getCode_Buffer(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Buffer>( new hardware::code::Buffer{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Spinors_staggered> getCode_Spinors_staggered(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Spinors_staggered>( new hardware::code::Spinors_staggered{*kernelParameters,deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Correlator_staggered> getCode_Correlator_staggered(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Correlator_staggered>( new hardware::code::Correlator_staggered{*kernelParameters, deviceIn} ) ;
		}
		std::unique_ptr<hardware::code::Fermions_staggered> getCode_Fermions_staggered(const hardware::Device * deviceIn) const
		{
			return std::unique_ptr<hardware::code::Fermions_staggered>( new hardware::code::Fermions_staggered{*kernelParameters, deviceIn} ) ;
		}
	private:
		const hardware::code::OpenClKernelParametersInterface * kernelParameters;
	};
}
