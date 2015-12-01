/*
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

#include "SpinorTester.hpp"
#include "GaugefieldTester.hpp"
#include "gaugefield.hpp"
#include "fermions.hpp"

//todo: is there a better solution to this?
//todo: can this be SpinorFillType instead of the vector? One only needs one anyway...
//todo: should this be called nonEvenOddFermionTestParameters?
struct FermionTestParameters : public NonEvenOddSpinorTestParameters, GaugefieldTestParameters
{
	FermionTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtentsIn, const SpinorFillTypes spinorFillTypesIn, const GaugefieldFillType gaugefieldFillTypesIn) :
		NonEvenOddSpinorTestParameters(referenceValuesIn, latticeExtentsIn, spinorFillTypesIn),
		GaugefieldTestParameters(referenceValuesIn, latticeExtentsIn, gaugefieldFillTypesIn) {};
	FermionTestParameters(const LatticeExtents lE, const SpinorFillTypes spinorFillTypesIn, const GaugefieldFillType gaugefieldFillTypesIn) :
		NonEvenOddSpinorTestParameters(defaultReferenceValues(), lE, spinorFillTypesIn),
		GaugefieldTestParameters(defaultReferenceValues(), lE, gaugefieldFillTypesIn) {};
	FermionTestParameters(const LatticeExtents lE, const SpinorFillTypes spinorFillTypesIn) :
		NonEvenOddSpinorTestParameters(defaultReferenceValues(), lE, spinorFillTypesIn),
		GaugefieldTestParameters(defaultReferenceValues(), lE, GaugefieldFillType::cold) {};
	FermionTestParameters() : NonEvenOddSpinorTestParameters(), GaugefieldTestParameters() {};
};

struct FermionTestParameters2 : public SpinorTestParameters, GaugefieldTestParameters
{
	FermionTestParameters2(const LatticeExtents lE, const SpinorFillTypes spinorFillTypesIn, const GaugefieldFillType gaugefieldFillTypesIn) :
		SpinorTestParameters(defaultReferenceValues(), lE, spinorFillTypesIn, false),
		GaugefieldTestParameters(defaultReferenceValues(), lE, gaugefieldFillTypesIn) {};
	FermionTestParameters2() : SpinorTestParameters(), GaugefieldTestParameters() {};
};

struct EvenOddFermionTestParameters : public EvenOddSpinorTestParameters, GaugefieldTestParameters
{
	EvenOddFermionTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn) :
		EvenOddSpinorTestParameters(referenceValuesIn, latticeExtentsIn, spinorFillTypeIn),
		GaugefieldTestParameters(referenceValuesIn, latticeExtentsIn, GaugefieldFillType::cold) {};
	EvenOddFermionTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtentsIn, const SpinorFillTypes spinorFillTypesIn,
			const GaugefieldFillType gaugefieldFillTypesIn = GaugefieldFillType::cold) :
			EvenOddSpinorTestParameters(referenceValuesIn, latticeExtentsIn, spinorFillTypesIn),
			GaugefieldTestParameters(referenceValuesIn, latticeExtentsIn, gaugefieldFillTypesIn) {};
};

struct FermionTester : public SpinorTester
{
	FermionTester(std::string kernelName, const ParameterCollection & pC, const FermionTestParameters & testParameters):
		SpinorTester(kernelName, pC, testParameters)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), device);
		const Matrixsu3 * gf_host = createGaugefield(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), testParameters.GaugefieldTestParameters::fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = SpinorTester::device->getFermionCode();
	}
	FermionTester(std::string kernelName, const ParameterCollection & pC, const EvenOddFermionTestParameters & testParameters):
		SpinorTester(kernelName, pC, testParameters)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), device);
		const Matrixsu3 * gf_host = createGaugefield(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), testParameters.GaugefieldTestParameters::fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = SpinorTester::device->getFermionCode();
	}

	FermionTester(std::string kernelName, const ParameterCollection & pC, const SpinorTestParameters & testParameters):
		SpinorTester(kernelName, pC, testParameters), gaugefieldBuffer(nullptr)
	{
		code = SpinorTester::device->getFermionCode();
	}
	
protected:
	const hardware::code::Fermions * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

//todo: this can probably be merged with FermionTester to be only one class (same for EvenOdd version)
struct FermionmatrixTester : public FermionTester
{
	FermionmatrixTester(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters) :
	FermionTester(kernelName, parameterCollection, testParameters)
	{
		in = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		out = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		in->load(createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));
	}
	~FermionmatrixTester()
	{
		calcSquarenormAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Plain<spinor> * in;
	const hardware::buffers::Plain<spinor> * out;
};

class FermionmatrixEvenOddTester : public FermionTester
{
public:
	FermionmatrixEvenOddTester(std::string kernelName, const ParameterCollection parameterCollection, const EvenOddFermionTestParameters testParameters) :
		FermionTester(kernelName, parameterCollection, testParameters)
	{
		in = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		out = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		in->load(createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));
	}
	~FermionmatrixEvenOddTester()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Spinor * in;
	const hardware::buffers::Spinor * out;
};

//@todo: rename nonTrivialMassParameter to nonTrivialParameter
const double nonTrivialMassParameter = 0.123456;

struct WilsonMassParameters
{
	WilsonMassParameters(const double kappaIn) : kappa(kappaIn){};
	const double kappa;
};

struct TwistedMassMassParameters : public WilsonMassParameters
{
	TwistedMassMassParameters(const double kappaIn, const double muIn):
		WilsonMassParameters(kappaIn), mu(muIn) {}

	const double mu;

	double getMubar() const { return 2.*kappa*mu; }
};

//todo: rename FermionTestParameters to NonEvenOddFermionTestParameters
template< class MassParameters>
 struct FermionMatrixTestParameters : public FermionTestParameters
{
	FermionMatrixTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
			const MassParameters massParametersIn, const ReferenceValues (* rV) (const int, const SpinorFillType, const GaugefieldFillType, const MassParameters)) :
		FermionTestParameters(rV( getSpinorfieldSize(latticeExtentsIn), spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn), latticeExtentsIn,
				SpinorFillTypes{spinorFillTypeIn}, gaugefieldFillTypeIn), massParameters(massParametersIn) {};
	FermionMatrixTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn,
			const MassParameters massParametersIn) :
		FermionTestParameters(defaultReferenceValues(), latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn}, gaugefieldFillTypeIn), massParameters(massParametersIn) {};
	FermionMatrixTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const MassParameters massParametersIn ) :
		FermionTestParameters(defaultReferenceValues(), latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn}, GaugefieldFillType::cold), massParameters(massParametersIn) {};
	const MassParameters massParameters;
};

struct EvenOddTwistedMassTestParameters: public EvenOddFermionTestParameters
{
	EvenOddTwistedMassTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn,
			const TwistedMassMassParameters massParametersIn, const ReferenceValues(*referenceValueCalculation) (int, TwistedMassMassParameters) ) :
		EvenOddFermionTestParameters(referenceValueCalculation( getEvenOddSpinorfieldSize(latticeExtentsIn), massParametersIn), latticeExtentsIn,
				SpinorFillTypes{spinorFillTypeIn}, GaugefieldFillType::cold), massParameters(massParametersIn) {};
	EvenOddTwistedMassTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gF,
			const TwistedMassMassParameters massParametersIn, const ReferenceValues(*referenceValueCalculation) (int, const GaugefieldFillType gF, TwistedMassMassParameters) ) :
		EvenOddFermionTestParameters(referenceValueCalculation( getEvenOddSpinorfieldSize(latticeExtentsIn), gF, massParametersIn), latticeExtentsIn,
				SpinorFillTypes{spinorFillTypeIn}, gF), massParameters(massParametersIn) {};

	const TwistedMassMassParameters massParameters;
};

struct FermionmatrixTester2a : public NonEvenOddSpinorTester
{
	FermionmatrixTester2a(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		NonEvenOddSpinorTester(kernelName, parameterCollection, testParameters, rV)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), device);
		const Matrixsu3 * gf_host = createGaugefield(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), testParameters.GaugefieldTestParameters::fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = SpinorTester::device->getFermionCode();

		in = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		out = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		in->load(createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));
	}
	~FermionmatrixTester2a()
	{
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Plain<spinor> * in;
	const hardware::buffers::Plain<spinor> * out;
	const hardware::code::Fermions * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct FermionmatrixTester2 : public FermionmatrixTester2a
{
	FermionmatrixTester2(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		FermionmatrixTester2a(kernelName, parameterCollection, testParameters, rV) {};
	~FermionmatrixTester2()
	{
		calcSquarenormAndStoreAsKernelResult(out);
	}
};

struct FermionmatrixEvenOddTester2 : public EvenOddSpinorTester
{
	FermionmatrixEvenOddTester2(std::string kernelName, const ParameterCollection parameterCollection, const FermionTestParameters testParameters, const ReferenceValues rV) :
		EvenOddSpinorTester(kernelName, parameterCollection, testParameters, rV)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), device);
		const Matrixsu3 * gf_host = createGaugefield(calculateGaugefieldSize(testParameters.SpinorTestParameters::latticeExtents), testParameters.GaugefieldTestParameters::fillType);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;

		code = SpinorTester::device->getFermionCode();

		in = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		out = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		in->load(createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));
	}
	~FermionmatrixEvenOddTester2()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Spinor * in;
	const hardware::buffers::Spinor * out;
	const hardware::code::Fermions * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};
