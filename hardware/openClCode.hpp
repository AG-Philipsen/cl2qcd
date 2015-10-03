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

#include "code/real.hpp"

#include "../meta/inputparameters.hpp"
#include <memory>

namespace hardware
{
	class OpenClCode
	{
	public:
		virtual ~OpenClCode(){};
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device *) const = 0;
	};

	class OpenClCode_fromMetaInputparameters final : public OpenClCode
	{
	public:
		OpenClCode_fromMetaInputparameters( const meta::Inputparameters & parametersIn ) : parameters(parametersIn) {};
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device * deviceIn) const override
		{
			std::unique_ptr<const hardware::code::Real> p1( new hardware::code::Real(parameters, deviceIn) ) ;
			return p1;
		}
	private:
		const meta::Inputparameters & parameters;
	};
}
