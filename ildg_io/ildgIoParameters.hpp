/** @file
 *
 * Interface for IldgIoParameters
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

#ifndef ILDGIOPARAMETERS_HPP_
#define ILDGIOPARAMETERS_HPP_

#include "../interfaceImplementations/latticesParameters.hpp"

class IldgIoParametersInterface
{
public:
	virtual ~IldgIoParametersInterface() {};
	virtual bool ignoreChecksumErrors() const = 0;
	virtual int getNumberOfElements() const = 0;
	virtual int getNs() const = 0;
	virtual int getNt() const = 0;
	virtual int getPrecision() const = 0;
	virtual double getKappa() const = 0;
	virtual double getMu() const = 0;
	virtual double getBeta() const = 0;
};

class Inputparameters : public IldgIoParametersInterface
{
public:
	Inputparameters (const physics::lattices::GaugefieldParametersInterface * parametersIn) : parameters(parametersIn) {} ;
	virtual bool ignoreChecksumErrors() const
	{
		return parameters->ignoreChecksumErrorsInIO();
	}
	virtual int getNumberOfElements() const
	{
		return parameters->getNumberOfElements();
	}
	virtual int getNs() const
	{
		return parameters->getNs();
	}
	virtual int getNt() const
	{
		return parameters->getNt();
	}
	virtual int getPrecision() const
	{
		return parameters->getPrecision();
	}
	virtual double getKappa() const
	{
		return parameters->getKappa();
	}
	virtual double getMu() const
	{
		return parameters->getMu();
	}
	virtual double getBeta() const
	{
		return parameters->getBeta();
	}
private:
	const physics::lattices::GaugefieldParametersInterface * parameters;
};

//TODO: this abstract class can be removed, the interface above is sufficient
class IldgIoParameters
{
public:
	IldgIoParameters(IldgIoParametersInterface * in) : itsClient( in ) {};
	virtual ~IldgIoParameters() {};
	virtual bool ignoreChecksumErrors() const = 0;
	virtual int getNumberOfElements() const = 0;
	virtual int getNt() const = 0;
	virtual int getNs() const = 0;
	virtual int getPrecision() const = 0;
	virtual double getKappa() const = 0;
	virtual double getBeta() const = 0;
	virtual double getMu() const = 0;
protected:
	//@Todo: Make this an instance??
	IldgIoParametersInterface * itsClient;
};

class IldgIoParameters_gaugefield: public IldgIoParameters
{
public:
	IldgIoParameters_gaugefield(IldgIoParametersInterface * in) : IldgIoParameters(in) {};
	virtual bool ignoreChecksumErrors() const
	{
		return itsClient->ignoreChecksumErrors();
	}
	virtual int getNumberOfElements() const
	{
		return itsClient->getNumberOfElements();
	}
	virtual int getNs() const
	{
		return itsClient->getNs();
	}
	virtual int getNt() const
	{
		return itsClient->getNt();
	}
	virtual int getPrecision() const
	{
		return itsClient->getPrecision();
	}
	virtual double getKappa() const
	{
		return itsClient->getKappa();
	}
	virtual double getMu() const
	{
		return itsClient->getMu();
	}
	virtual double getBeta() const
	{
		return itsClient->getBeta();
	}
};

IldgIoParameters_gaugefield createIldgIoParameters(const meta::Inputparameters * parametersIn);

#endif /* ILDGIOPARAMETERS_HPP_ */
