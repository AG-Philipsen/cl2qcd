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

#include "../../common_header_files/types.h"
#include <iostream>

class LatticeObjectParametersInterface {
public:
	virtual ~LatticeObjectParametersInterface () {};
	virtual int getNs() const = 0;
	virtual int getNt() const = 0;
	virtual int getPrecision() const = 0;
	virtual bool ignoreChecksumErrorsInIO() const = 0;
	virtual int getNumberOfElements() const = 0;
	virtual double getKappa() const = 0;
	virtual double getMu() const = 0;
	virtual double getBeta() const = 0;
	virtual common::startcondition getStartcondition() const = 0;
	virtual std::string getNamePrefix() const = 0;
	virtual std::string getNamePostfix() const = 0;
	virtual int getNumberOfDigitsInName() const = 0;
	virtual int getSmearingSteps() const = 0;
	virtual std::string getSourcefileName() const = 0;
};

#include "../../meta/inputparameters.hpp"
#include "../../meta/util.hpp"

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
	virtual bool ignoreChecksumErrorsInIO() const
	{
		return parameters->get_ignore_checksum_errors();
	}
	virtual int getNumberOfElements() const
	{
		return meta::get_vol4d(*parameters) * NDIM;
	}
	virtual double getKappa() const
	{
		return parameters->get_kappa();
	}
	virtual double getMu() const
	{
		return parameters->get_mu();
	}
	virtual double getBeta() const
	{
		return parameters->get_beta();
	}
	virtual common::startcondition getStartcondition() const
	{
		return parameters->get_startcondition();
	}
	virtual std::string getNamePrefix() const
	{
		return parameters->get_config_prefix();
	}
	virtual std::string getNamePostfix() const
	{
		return parameters->get_config_postfix();
	}
	virtual int getNumberOfDigitsInName() const
	{
		return parameters->get_config_number_digits();
	}
	virtual int getSmearingSteps() const
	{
		return parameters->get_rho_iter();
	}
	virtual std::string getSourcefileName() const
	{
		return parameters->get_sourcefile();
	}

private:
	const meta::Inputparameters * parameters;
};

class SpinorfieldParametersImplementation final {
    public:
        //TODO: uncomment following line
        //SpinorfieldParametersImplementation() = delete;
        SpinorfieldParametersImplementation(const meta::Inputparameters& paramsIn ): parameters(paramsIn) {};
        int getNt() const
        {
            return parameters.get_ntime();
        }
        int getNs() const
        {
            return parameters.get_nspace();
        }

    private:
       const meta::Inputparameters& parameters;


};



