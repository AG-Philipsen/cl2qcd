/*
 * Copyright 2012, 2013, 2014 Christopher Pinke, Matthias Bach
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_FERMIONS

#include "gaugefield.hpp"
#include "fermions.hpp"
#include "../../physics/lattices/gaugefield.hpp"

#include "FermionTester.hpp"

/** todo:
 * - the way the reference values are passed is not nice, it should be somehow done automatically.
 * 		Especially that one has to call the lattice volume explicitely is not nice
 * - Using templates makes the actual test commands lenghty. Perhaps it is a good idea to create an
 * 	in-between fct. which passes the arguments to the template.
 */

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		FermionTester tester("build", "build_input_1");
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
		FermionTester tester("build", "build_input_2");
	}

	BOOST_AUTO_TEST_CASE( BUILDFROMPARAMETERS )
	{
		const hardware::HardwareParametersMockup hardwareParameters(4,4);
		const hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(4,4);
		ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
		const FermionTestParameters testParameters;
		BOOST_CHECK_NO_THROW( FermionTester( "build all kernels", parameterCollection, testParameters) );
	}

BOOST_AUTO_TEST_SUITE_END()

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
	FermionmatrixEvenOddTester(std::string kernelName, std::string inputfile) :
	FermionTester(kernelName, inputfile, 1)
	{
		in = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		out = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		in->load(createSpinorfield(spinorfieldEvenOddElements));
		out->load(createSpinorfield(spinorfieldEvenOddElements));
	}
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

const ReferenceValues calculateReferenceValues_mWilson(const int latticeVolume, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn)
{
	if( massParametersIn.kappa == 0. and spinorFillTypeIn == SpinorFillType::ascendingComplex)
	{
		return ReferenceValues{latticeVolume * sumOfIntegersSquared(24)};
	}
	else if (massParametersIn.kappa == nonTrivialMassParameter)
	{
		if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::cold )
		{
			return ReferenceValues{latticeVolume * 0.7476023296000095 };
		}
		else if (spinorFillTypeIn == SpinorFillType::one and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
		{
			return ReferenceValues{latticeVolume * 7.03532241159178};
		}
		else if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial)
		{
			return ReferenceValues{latticeVolume * 3508.912843225604};
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmMinus(const int latticeVolume, const SpinorFillType spinorFillTypeIn,
		const GaugefieldFillType gaugefieldFillTypeIn, const TwistedMassMassParameters massParametersIn)
{
	if( massParametersIn.kappa == nonTrivialMassParameter and spinorFillTypeIn == SpinorFillType::ascendingComplex )
	{
		if( massParametersIn.mu == 0. )
		{
			return calculateReferenceValues_mWilson(latticeVolume, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
		}
		else if (massParametersIn.mu == nonTrivialMassParameter and gaugefieldFillTypeIn == GaugefieldFillType::cold )
		{
			return ReferenceValues{latticeVolume * 5.300678101577363};
		}
		else if (massParametersIn.mu == nonTrivialMassParameter and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
		{
			return ReferenceValues{latticeVolume * 3513.465918997579};
		}
	}
	return defaultReferenceValues();
}

//todo: these two can be combined with the reference fcts. as template arguments
struct MWilsonTestParameters : public FermionTestParameters
{
	MWilsonTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn) :
		FermionTestParameters(calculateReferenceValues_mWilson( getSpinorfieldSize(latticeExtentsIn), spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn), latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn}, gaugefieldFillTypeIn), massParameters(massParametersIn) {};
	const WilsonMassParameters massParameters;
};

struct TwistedMassTestParameters : public FermionTestParameters
{
	TwistedMassTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const TwistedMassMassParameters massParametersIn) :
		FermionTestParameters(calculateReferenceValues_mTmMinus( getSpinorfieldSize(latticeExtentsIn), spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn), latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn}, gaugefieldFillTypeIn), massParameters(massParametersIn) {};
	const TwistedMassMassParameters massParameters;
};

template<typename TesterClass, typename MassParameters, typename KernelParameterMockup, typename TestParameters>
void callTest( const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const MassParameters massParametersIn)
{
	TestParameters parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	KernelParameterMockup kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE( M_WILSON )

	class MWilsonTester : public FermionmatrixTester
	{
	public:
		MWilsonTester(const ParameterCollection parameterCollection, const MWilsonTestParameters & testParameters) :
		FermionmatrixTester("m_wilson", parameterCollection, testParameters)
		{
			code->M_wilson_device(in, out,  gaugefieldBuffer, testParameters.massParameters.kappa );
		}
	};

	BOOST_AUTO_TEST_CASE( M_WILSON_1)
	{
		callTest<MWilsonTester, WilsonMassParameters,hardware::code::OpenClKernelParametersMockupForSpinorTests,MWilsonTestParameters>(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{0.});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_2)
	{
		callTest<MWilsonTester, WilsonMassParameters,hardware::code::OpenClKernelParametersMockupForSpinorTests,MWilsonTestParameters>(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialMassParameter});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_3)
	{
		callTest<MWilsonTester, WilsonMassParameters,hardware::code::OpenClKernelParametersMockupForSpinorTests,MWilsonTestParameters>(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialMassParameter});
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_4)
	{
		callTest<MWilsonTester, WilsonMassParameters,hardware::code::OpenClKernelParametersMockupForSpinorTests,MWilsonTestParameters>(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  )

	class MTmMinusTester : public FermionmatrixTester
	{
	public:
		MTmMinusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & testParameters) :
			FermionmatrixTester("m_tm_minus", parameterCollection, testParameters)
		{
			code->M_tm_minus_device(in, out,  gaugefieldBuffer, testParameters.massParameters.kappa, testParameters.massParameters.getMubar() );
		}
	};

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 )
	{
		callTest<MTmMinusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns12, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialMassParameter, 0.});
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_2 )
	{
		callTest<MTmMinusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns4, nt16}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialMassParameter, nonTrivialMassParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_3 )
	{
		callTest<MTmMinusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{nonTrivialMassParameter, nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

//todo: add own parameters
//todo: test should not pass using M_tm_minus ref. values!
BOOST_AUTO_TEST_SUITE( M_TM_PLUS )

	class MTmPlusTester : public FermionmatrixTester
	{
	public:
		MTmPlusTester(const ParameterCollection parameterCollection, const TwistedMassTestParameters & testParameters) :
			FermionmatrixTester("m_tm_plus", parameterCollection, testParameters)
		{
			code->M_tm_plus_device(in, out,  gaugefieldBuffer, testParameters.massParameters.kappa, testParameters.massParameters.getMubar() );
		}
	};

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_1 )
	{
		callTest<MTmPlusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns12, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialMassParameter, 0.});
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_2 )
	{
		callTest<MTmPlusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns4, nt16}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{nonTrivialMassParameter, nonTrivialMassParameter});
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_3 )
	{
		callTest<MTmPlusTester, TwistedMassMassParameters,hardware::code::OpenClKernelParametersMockupForTwistedMass,TwistedMassTestParameters>(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{nonTrivialMassParameter, nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 )

	ReferenceValues calculateReferenceValues_gamma5(const int latticeVolume, const SpinorFillType spinorFillTypeIn)
	{
		return (spinorFillTypeIn == SpinorFillType::ascendingComplex) ? ReferenceValues{latticeVolume*( sumOfIntegers(1, 12, 1) - sumOfIntegers(13, 24, 1) )} : defaultReferenceValues() ;
	}

	//todo: merge these two
	struct Gamma5Tester : public FermionTester
	{
		Gamma5Tester(const ParameterCollection parameterCollection, const NonEvenOddSpinorTestParameters & testParameters) :
			FermionTester("gamma5", parameterCollection, testParameters)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			spinor * sf_in;
			sf_in = new spinor[spinorfieldElements];

			in.load( createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)) );
			code->gamma5_device(&in);
			in.dump(sf_in);
			kernelResult.at(0) = count_sf(sf_in, spinorfieldElements);

			delete sf_in;
		}
	};

	struct Gamma5EvenOddTester : public FermionTester
	{
		Gamma5EvenOddTester(const ParameterCollection parameterCollection, const EvenOddSpinorTestParameters & testParameters) :
			FermionTester("gamma5_eo", parameterCollection, testParameters)
		{
			const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
			spinor * sf_in;
			sf_in = new spinor[spinorfieldEvenOddElements];

			in.load( createSpinorfield(testParameters.SpinorTestParameters::fillTypes.at(0)) );
			code->gamma5_eo_device(&in);
			in.dump(sf_in);
			kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);

			delete sf_in;
		}
	};

	template<class TesterClass, class ParameterClass>
	void runTest(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn)
	{
		ParameterClass parametersForThisTest(referenceValuesIn, latticeExtentsIn, spinorFillTypeIn);
		hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
		hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
		ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
		TesterClass tester(parameterCollection, parametersForThisTest);
	}

	BOOST_AUTO_TEST_CASE( GAMMA5_NONEO )
	{
		runTest<Gamma5Tester, NonEvenOddSpinorTestParameters>(calculateReferenceValues_gamma5(calculateSpinorfieldSize( LatticeExtents{ns4,nt4}), SpinorFillType::ascendingComplex), LatticeExtents{ns4,nt4}, SpinorFillType::ascendingComplex);
	}

	BOOST_AUTO_TEST_CASE( GAMMA5_EO )
	{
		runTest<Gamma5EvenOddTester, EvenOddSpinorTestParameters>(calculateReferenceValues_gamma5(calculateEvenOddSpinorfieldSize( LatticeExtents{ns4,nt4}), SpinorFillType::ascendingComplex), LatticeExtents{ns4,nt4}, SpinorFillType::ascendingComplex);
	}

BOOST_AUTO_TEST_SUITE_END()

template<class TesterClass, class ParameterClass>
void runTest(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParameters)
{
	ParameterClass parametersForThisTest(latticeExtentsIn, spinorFillTypeIn, massParameters);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

ReferenceValues calculateReferenceValues_mTmSitediagonal(const int latticeVolume, const TwistedMassMassParameters massParametersIn)
{
	if (massParametersIn.kappa == nonTrivialMassParameter and massParametersIn.mu == nonTrivialMassParameter)
	{
		return ReferenceValues{ latticeVolume * 4904.553075771979};
	}
	return defaultReferenceValues();
}

//todo: remove ARG_DEF from all the tm diagonal kernel fcts.!
BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL )

	struct MTmSitediagonalTestParameters: public EvenOddFermionTestParameters
	{
		MTmSitediagonalTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillType spinorFillTypeIn, const TwistedMassMassParameters massParametersIn) :
			EvenOddFermionTestParameters(calculateReferenceValues_mTmSitediagonal( getEvenOddSpinorfieldSize(latticeExtentsIn), massParametersIn), latticeExtentsIn, SpinorFillTypes{spinorFillTypeIn}, GaugefieldFillType::cold), massParameters(massParametersIn) {};
		const TwistedMassMassParameters massParameters;
	};

	struct MTmSitediagonalTester: public FermionmatrixEvenOddTester
	{
		MTmSitediagonalTester(const ParameterCollection parameterCollection, const MTmSitediagonalTestParameters & testParameters):
			FermionmatrixEvenOddTester("m_tm_sitediagonal", parameterCollection, testParameters)
			{
				code->M_tm_sitediagonal_device( in, out, testParameters.massParameters.getMubar());
			}
	};

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_1)
	{
		runTest<MTmSitediagonalTester, MTmSitediagonalTestParameters>(LatticeExtents{ns4,nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{nonTrivialMassParameter, nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL )

	class MTmInverseSitediagonalTester: public FermionmatrixEvenOddTester
	{
	public:
		MTmInverseSitediagonalTester(std::string inputfile):
			FermionmatrixEvenOddTester("m_tm_inverse_sitediagonal", inputfile)
			{
				code->M_tm_inverse_sitediagonal_device( in, out);
			}
	};

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1)
	{
		MTmInverseSitediagonalTester tester("m_tm_inverse_sitediagonal_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_2)
	{
		MTmInverseSitediagonalTester tester("m_tm_inverse_sitediagonal_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_3)
	{
		MTmInverseSitediagonalTester tester("m_tm_inverse_sitediagonal_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS )

	class MTmSitediagonalMinusTester: public FermionmatrixEvenOddTester
	{
	public:
		MTmSitediagonalMinusTester(std::string inputfile):
			FermionmatrixEvenOddTester("m_tm_sitediagonal_minus", inputfile)
			{
				code->M_tm_sitediagonal_minus_device( in, out);
			}
	};

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1)
	{
		MTmSitediagonalMinusTester tester("m_tm_sitediagonal_minus_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_2)
	{
		MTmSitediagonalMinusTester tester("m_tm_sitediagonal_minus_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_3)
	{
		MTmSitediagonalMinusTester tester("m_tm_sitediagonal_minus_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL_MINUS)

	class MTmInverseSitediagonalMinusTester: public FermionmatrixEvenOddTester
	{
	public:
		MTmInverseSitediagonalMinusTester(std::string inputfile):
			FermionmatrixEvenOddTester("m_tm_inverse_sitediagonal_minus", inputfile)
			{
				code->M_tm_inverse_sitediagonal_minus_device( in, out);
			}
	};

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1)
	{
		MTmInverseSitediagonalMinusTester tester("m_tm_inverse_sitediagonal_minus_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_2)
	{
		MTmInverseSitediagonalMinusTester tester("m_tm_inverse_sitediagonal_minus_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_3)
	{
		MTmInverseSitediagonalMinusTester tester("m_tm_inverse_sitediagonal_minus_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_EO )

	class DslashEvenOddTester: public FermionmatrixEvenOddTester
	{
	public:
		DslashEvenOddTester(std::string inputfile):
			FermionmatrixEvenOddTester("dslash_eo", inputfile)
			{
				evenOrOdd ? 
					code->dslash_eo_device( in, out, this->getGaugefieldBuffer(), EVEN, parameters->get_kappa() ) :
					code->dslash_eo_device( in, out, this->getGaugefieldBuffer(), ODD, parameters->get_kappa() );
			}
	};

	BOOST_AUTO_TEST_CASE( DSLASH_EO_1)
	{
		DslashEvenOddTester tester("dslash_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_2)
	{
		DslashEvenOddTester tester("dslash_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_3)
	{
		DslashEvenOddTester tester("dslash_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_4)
	{
		DslashEvenOddTester tester("dslash_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_5)
	{
		DslashEvenOddTester tester("dslash_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_6)
	{
		DslashEvenOddTester tester("dslash_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_7)
	{
		DslashEvenOddTester tester("dslash_eo_input_7");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_8)
	{
		DslashEvenOddTester tester("dslash_eo_input_8");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_9)
	{
		DslashEvenOddTester tester("dslash_eo_input_9");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_10)
	{
		DslashEvenOddTester tester("dslash_eo_input_10");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_11)
	{
		DslashEvenOddTester tester("dslash_eo_input_11");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_12)
	{
		DslashEvenOddTester tester("dslash_eo_input_12");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_13)
	{
		DslashEvenOddTester tester("dslash_eo_input_13");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_14)
	{
		DslashEvenOddTester tester("dslash_eo_input_14");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_15)
	{
		DslashEvenOddTester tester("dslash_eo_input_15");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_16)
	{
		DslashEvenOddTester tester("dslash_eo_input_16");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_17)
	{
		DslashEvenOddTester tester("dslash_eo_input_17");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_18)
	{
		DslashEvenOddTester tester("dslash_eo_input_18");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_19)
	{
		DslashEvenOddTester tester("dslash_eo_input_19");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_20)
	{
		DslashEvenOddTester tester("dslash_eo_input_20");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_21)
	{
		DslashEvenOddTester tester("dslash_eo_input_21");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_22)
	{
		DslashEvenOddTester tester("dslash_eo_input_22");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_23)
	{
		DslashEvenOddTester tester("dslash_eo_input_23");
	}

	BOOST_AUTO_TEST_CASE( DSLASH_EO_24)
	{
		DslashEvenOddTester tester("dslash_eo_input_24");
	}

BOOST_AUTO_TEST_SUITE_END()

class MFermionEvenOddComparator : public FermionTester
{
public:
	MFermionEvenOddComparator(std::string kernelName, std::string inputfile):
		FermionTester(kernelName + " in EvenOdd and non-EvenOdd formulation", inputfile, 2)
	{
		in_eo1 = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		in_eo2 = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		out_eo = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		in_noneo = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		out_noneo = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		intermediate1 = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		intermediate2 = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		intermediate3 = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		
		in_eo1->load( createSpinorfield(spinorfieldEvenOddElements, 123456 ));
		in_eo2->load( createSpinorfield(spinorfieldEvenOddElements, 78910 ));
		in_noneo->load ( createSpinorfield(spinorfieldElements, 123456) );

		if(useRandom)
		{
			evenOrOdd ? 
				SpinorTester::code->convert_from_eoprec_device(in_eo1, in_eo2, in_noneo) :
				SpinorTester::code->convert_to_eoprec_device(in_eo1, in_eo2, in_noneo);
		}
		
		hmc_complex minusone_tmp = { -1., 0.};
		minusone = new hardware::buffers::Plain<hmc_complex>(1, device);
		minusone->load(&minusone_tmp);
		
		referenceValue[1] = referenceValue[0];
	}
	~MFermionEvenOddComparator()
	{
		delete in_eo1;
		delete in_eo2;
		delete out_eo;
		delete in_noneo;
		delete out_noneo,
		delete intermediate1;
		delete intermediate2;
		delete intermediate3;
		
		kernelResult[0] = resultEvenOdd;
		kernelResult[1] = resultNonEvenOdd;
	}
	
protected:
		double resultEvenOdd;
		double resultNonEvenOdd;
	
		const hardware::buffers::Spinor * in_eo1;
		const hardware::buffers::Spinor * in_eo2;
		const hardware::buffers::Plain<spinor> * out_eo;
		const hardware::buffers::Plain<spinor> * in_noneo;
		const hardware::buffers::Plain<spinor> * out_noneo;
		const hardware::buffers::Spinor * intermediate1;
		const hardware::buffers::Spinor * intermediate2;
		const hardware::buffers::Spinor * intermediate3;
		hardware::buffers::Plain<hmc_complex> * minusone;
};




BOOST_AUTO_TEST_SUITE(M_WILSON_COMPARE_NONEO_EO )

	class MWilsonEvenOddComparator : public MFermionEvenOddComparator
	{
	public:
		MWilsonEvenOddComparator(std::string inputfile):
			MFermionEvenOddComparator("m_wilson", inputfile)
			{
				code->M_wilson_device(in_noneo, out_noneo,  this->getGaugefieldBuffer(), parameters->get_kappa());
				SpinorTester::code->set_float_to_global_squarenorm_device(out_noneo, doubleBuffer);
				doubleBuffer->dump(&resultNonEvenOdd);
				
				//suppose in1 is the even, in2 the odd input vector
				//now calc out_tmp_eo1 = (1 in1 + D_eo in2)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate1);
				code->dslash_eo_device(in_eo2, intermediate1, this->getGaugefieldBuffer(), EO, parameters->get_kappa());
				SpinorTester::code->saxpy_eoprec_device(intermediate1, in_eo1, minusone, intermediate1);

				//now calc out_tmp_eo2 = ( 1 in2 + D_oe in1)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate2);
				code->dslash_eo_device(in_eo1, intermediate2, this->getGaugefieldBuffer(), OE, parameters->get_kappa());
				SpinorTester::code->saxpy_eoprec_device(intermediate2, in_eo2, minusone, intermediate2);

				SpinorTester::code->convert_from_eoprec_device(intermediate1, intermediate2, out_eo);
				
				SpinorTester::code->set_float_to_global_squarenorm_device(out_eo, doubleBuffer);
				doubleBuffer->dump(&resultEvenOdd);
			}
	};

	BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_1)
	{
		MWilsonEvenOddComparator tester("m_wilson_compare_noneo_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_2)
	{
		MWilsonEvenOddComparator tester("m_wilson_compare_noneo_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_3)
	{
		MWilsonEvenOddComparator tester("m_wilson_compare_noneo_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_4)
	{
		MWilsonEvenOddComparator tester("m_wilson_compare_noneo_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_5)
	{
		MWilsonEvenOddComparator tester("m_wilson_compare_noneo_eo_input_5");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_COMPARE_NONEO_EO )

	class MTmEvenOddComparator : public MFermionEvenOddComparator
	{
	public:
		MTmEvenOddComparator(std::string inputfile):
			MFermionEvenOddComparator("m_tm_plus", inputfile)
			{
				code->M_tm_plus_device(in_noneo, out_noneo,  this->getGaugefieldBuffer(), parameters->get_kappa(), meta::get_mubar(*parameters));
				SpinorTester::code->set_float_to_global_squarenorm_device(out_noneo, doubleBuffer);
				doubleBuffer->dump(&resultNonEvenOdd);
				
				//suppose in1 is the even, in2 the odd input vector
				//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate1);
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate2);

				code->dslash_eo_device(in_eo2, intermediate1, this->getGaugefieldBuffer(), EO, parameters->get_kappa());
				code->M_tm_sitediagonal_device(in_eo1, intermediate2,  meta::get_mubar(*parameters));

				SpinorTester::code->saxpy_eoprec_device(intermediate1, intermediate2, minusone, intermediate1);

				//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate2);
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate3);

				code->dslash_eo_device(in_eo1, intermediate2, this->getGaugefieldBuffer(), OE, parameters->get_kappa());
				code->M_tm_sitediagonal_device(in_eo2, intermediate3,  meta::get_mubar(*parameters));

				SpinorTester::code->saxpy_eoprec_device(intermediate2, intermediate3, minusone, intermediate2);

				//now, both output vectors have to be converted back to noneo
				SpinorTester::code->convert_from_eoprec_device(intermediate1, intermediate2, out_eo);
				SpinorTester::code->set_float_to_global_squarenorm_device(out_eo, doubleBuffer);
				doubleBuffer->dump(&resultEvenOdd);
			}
	};

	BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_1)
	{
		MTmEvenOddComparator tester("m_tm_compare_noneo_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_2)
	{
		MTmEvenOddComparator tester("m_tm_compare_noneo_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_3)
	{
		MTmEvenOddComparator tester("m_tm_compare_noneo_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_4)
	{
		MTmEvenOddComparator tester("m_tm_compare_noneo_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_5)
	{
		MTmEvenOddComparator tester("m_tm_compare_noneo_eo_input_5");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_MINUS_COMPARE_NONEO_EO )

	class MTmMinusEvenOddComparator : public MFermionEvenOddComparator
	{
	public:
		MTmMinusEvenOddComparator(std::string inputfile):
			MFermionEvenOddComparator("m_tm_minus", inputfile)
			{
				code->M_tm_minus_device(in_noneo, out_noneo,  this->getGaugefieldBuffer(), parameters->get_kappa(), meta::get_mubar(*parameters));
				SpinorTester::code->set_float_to_global_squarenorm_device(out_noneo, doubleBuffer);
				doubleBuffer->dump(&resultNonEvenOdd);
				
				//suppose in1 is the even, in2 the odd input vector
				//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate1);
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate2);

				code->dslash_eo_device(in_eo2, intermediate1, this->getGaugefieldBuffer(), EO, parameters->get_kappa());
				code->M_tm_sitediagonal_minus_device(in_eo1, intermediate2,  meta::get_mubar(*parameters));

				SpinorTester::code->saxpy_eoprec_device(intermediate1, intermediate2, minusone, intermediate1);

				//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate2);
				SpinorTester::code->set_zero_spinorfield_eoprec_device(intermediate3);

				code->dslash_eo_device(in_eo1, intermediate2, this->getGaugefieldBuffer(), OE, parameters->get_kappa());
				code->M_tm_sitediagonal_minus_device(in_eo2, intermediate3,  meta::get_mubar(*parameters));

				SpinorTester::code->saxpy_eoprec_device(intermediate2, intermediate3, minusone, intermediate2);

				//now, both output vectors have to be converted back to noneo
				SpinorTester::code->convert_from_eoprec_device(intermediate1, intermediate2, out_eo);
				SpinorTester::code->set_float_to_global_squarenorm_device(out_eo, doubleBuffer);
				doubleBuffer->dump(&resultEvenOdd);
			}
	};

	BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_1)
	{
		MTmMinusEvenOddComparator tester("m_tm_minus_compare_noneo_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_2)
	{
		MTmMinusEvenOddComparator tester("m_tm_minus_compare_noneo_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_3)
	{
		MTmMinusEvenOddComparator tester("m_tm_minus_compare_noneo_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_4)
	{
		MTmMinusEvenOddComparator tester("m_tm_minus_compare_noneo_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_5)
	{
		MTmMinusEvenOddComparator tester("m_tm_minus_compare_noneo_eo_input_5");
	}

BOOST_AUTO_TEST_SUITE_END()

