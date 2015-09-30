/** @file
 * Declaration of the physics::lattices::Gaugefield class
 *
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

#pragma once

class LatticeObjectParametersInterface {
public:
	virtual ~LatticeObjectParametersInterface () {};
	virtual int getNs() const = 0;
	virtual int getNt() const = 0;
	virtual int getPrecision() const = 0;
};

#include "../../meta/inputparameters.hpp"

#include <iostream>
class LatticeObjectParametersImplementation : public LatticeObjectParametersInterface{
public:
	LatticeObjectParametersImplementation() = delete;
	LatticeObjectParametersImplementation(const meta::Inputparameters * paramsIn ): parameters(paramsIn) {};
	virtual ~LatticeObjectParametersImplementation () {};
	virtual int getNs() const
	{
		return parameters->get_nspace();
	}
	virtual int getNt() const
	{
		return parameters->get_ntime();
	}
	virtual int getPrecision() const
	{
		return parameters->get_precision();
	}
private:
	const meta::Inputparameters * parameters;
};
