/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_SPINORS_STAGGERED

#include "SpinorStaggeredTester.hpp"
#include "Kolmogorov_Smirnov.h"
#include "Normal_RNG_tests.h"

#include "../../host_functionality/logger.hpp"

const ReferenceValues calculateReferenceValues_squarenorm(const int latticeVolume, const SpinorFillTypes fillTypesIn)
{
	switch( fillTypesIn.at(0) )
	{
		case SpinorFillType::one :
		{
			return ReferenceValues{latticeVolume * 3.};
		}
		case SpinorFillType::ascendingComplex:
		{
			return ReferenceValues{latticeVolume * sumOfIntegersSquared(6)};
		}
		default:
		{
			return defaultReferenceValues();
		}
	}
}

const ReferenceValues calculateReferenceValues_scalarProduct(const int latticeVolume, const SpinorFillTypes fillTypesIn)
{
	if(fillTypesIn.at(0) == SpinorFillType::one and fillTypesIn.at(1) == SpinorFillType::one)
	{
		return ReferenceValues{latticeVolume * 3., 0.};
	}
	else if( fillTypesIn.at(0) == SpinorFillType::one and fillTypesIn.at(1) == SpinorFillType::ascendingComplex )
	{
		return ReferenceValues{latticeVolume * sumOfIntegers(1,5,2), latticeVolume * sumOfIntegers(2,6,2)};
	}
	else if( fillTypesIn.at(0) == SpinorFillType::ascendingComplex and fillTypesIn.at(1) == SpinorFillType::one )
	{
		return ReferenceValues{latticeVolume * sumOfIntegers(1,5,2), -latticeVolume * sumOfIntegers(2,6,2)};
	}
	else if ( fillTypesIn.at(0) == SpinorFillType::ascendingComplex and fillTypesIn.at(1) == SpinorFillType::ascendingComplex  )
	{
		return ReferenceValues{latticeVolume * sumOfIntegersSquared(6), 0.};
	}
	else
		return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_cold(const bool isEvenOdd)
{
	return (isEvenOdd) ? ReferenceValues{0.5} : ReferenceValues{1.};
}

const ReferenceValues calculateReferenceValues_zero()
{
	return ReferenceValues{0.};
}

const ReferenceValues calculateReferenceValues_sax (const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{ (alphaIn.at(0).im * alphaIn.at(0).im + alphaIn.at(0).re * alphaIn.at(0).re) * latticeVolume * sumOfIntegersSquared(6)};
}

const ReferenceValues calculateReferenceValues_saxpy(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{1. + alphaIn.at(0).re, 0. + alphaIn.at(0).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_saxpbypz(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{1. + alphaIn.at(0).re + alphaIn.at(1).re, 0. + alphaIn.at(0).im + alphaIn.at(1).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_convert_eo(const int latticeVolume, const bool fillEvenSites)
{
	double nonTrivialValue = latticeVolume * 3.;
	return ReferenceValues {fillEvenSites ? nonTrivialValue : 0, fillEvenSites ? 0. : nonTrivialValue};
}

const ReferenceValues calculateReferenceValues_convertFromEvenOdd(const int latticeVolume)
{
	return ReferenceValues {latticeVolume * 12. / 2.};
}

struct LinearCombinationTestParameters
{
	LinearCombinationTestParameters(const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn) : coefficients(coefficientsIn), numberOfSpinors(numberOfSpinorsIn) {};
	LinearCombinationTestParameters(const size_t numberOfSpinorsIn) : coefficients(ComplexNumbers{}), numberOfSpinors(numberOfSpinorsIn) {};
	LinearCombinationTestParameters() : coefficients(ComplexNumbers{{1.,0.}}), numberOfSpinors(1) {};
	const ComplexNumbers coefficients;
	const NumberOfSpinors numberOfSpinors;
};

struct NonEvenOddLinearCombinationTestParameters : public NonEvenOddSpinorStaggeredTestParameters, LinearCombinationTestParameters
{
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn) :
		NonEvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn), LinearCombinationTestParameters{coefficientsIn, numberOfSpinorsIn} {};
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const size_t numberOfSpinorsIn, const ComparisonType typeOfComparisonIn):
		NonEvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisonIn), LinearCombinationTestParameters{numberOfSpinorsIn}{};
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		NonEvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn) {};
};

struct EvenOddLinearCombinationTestParameters : public EvenOddSpinorStaggeredTestParameters, LinearCombinationTestParameters
{
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn) :
		EvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn), LinearCombinationTestParameters{coefficientsIn, numberOfSpinorsIn} {};
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const size_t numberOfSpinorsIn, const ComparisonType typeOfComparisonIn):
		EvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisonIn), LinearCombinationTestParameters{numberOfSpinorsIn}{};
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		EvenOddSpinorStaggeredTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn) {};
};

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass, typename AdditionalArgument> void performTest(LatticeExtents latticeExtendsIn, const AdditionalArgument addArg )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, addArg);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, fillTypesIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(LatticeExtents latticeExtendsIn, const ComplexNumbers alphaIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, SpinorFillTypes{SpinorFillType::ascendingComplex}, alphaIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(LatticeExtents latticeExtendsIn, const bool fillEvenSitesIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, fillEvenSitesIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(LatticeExtents latticeExtendsIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisonIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, typeOfComparisonIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	   BOOST_CHECK_NO_THROW(SpinorStaggeredTester spinorStaggeredTester("build all kernels",
									     "spinors_staggered_build_input_1", 0));
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	   BOOST_CHECK_NO_THROW(SpinorStaggeredTester spinorStaggeredTester("build all kernels", 
									      "spinors_staggered_build_input_2", 0));
	}

	BOOST_AUTO_TEST_CASE( BUILDFROMPARAMETERS )
	{
		const hardware::HardwareParametersMockup hardwareParameters(4,4);
		const hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(4,4);
		ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
		const SpinorStaggeredTestParameters testParameters;
		BOOST_CHECK_NO_THROW( SpinorStaggeredTester( "build all kernels", parameterCollection, testParameters) );
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

class LinearCombinationTester: public SpinorStaggeredTester
{
public:
	LinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const NonEvenOddLinearCombinationTestParameters testParameters):
		SpinorStaggeredTester(kernelName, parameterCollection, testParameters)
	{
		loadCoefficients(testParameters);
	}
	LinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const EvenOddLinearCombinationTestParameters testParameters):
		SpinorStaggeredTester(kernelName, parameterCollection, testParameters)
	{
		loadCoefficients(testParameters);
	}
protected:
	std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
private:
	void loadCoefficients(const LinearCombinationTestParameters & testParameters)
	{
		for (auto coefficient : testParameters.coefficients)
		{
			complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
			complexNums.back()->load(&coefficient);
		}
	}
};

class NonEvenOddLinearCombinationTester: public LinearCombinationTester
{
public:
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const NonEvenOddLinearCombinationTestParameters testParameters):
		LinearCombinationTester(kernelName, parameterCollection, testParameters)
	{
		for( size_t number = 0; number < testParameters.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::Plain<su3vec>(spinorfieldElements, device));
			(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
		}
	}
protected:
	std::vector<const hardware::buffers::Plain<su3vec> *> spinorfields;
	const hardware::buffers::Plain<su3vec> * getOutSpinor() const
	{
		return spinorfields.back();
	}
};

struct NonEvenOddLinearCombinationTesterWithSquarenormAsResult : public NonEvenOddLinearCombinationTester
{
	NonEvenOddLinearCombinationTesterWithSquarenormAsResult(const std::string kernelName, const ParameterCollection parameterCollection, const NonEvenOddLinearCombinationTestParameters testParameters) :
		NonEvenOddLinearCombinationTester(kernelName, parameterCollection, testParameters) {};

	~NonEvenOddLinearCombinationTesterWithSquarenormAsResult()
	{
		calcSquarenormAndStoreAsKernelResult(getOutSpinor());
	}
};

class EvenOddLinearCombinationTester: public LinearCombinationTester
{
public:
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection & parameterCollection, const EvenOddLinearCombinationTestParameters & testParameters):
		LinearCombinationTester(kernelName, parameterCollection, testParameters)
	{
		for( size_t number = 0; number < testParameters.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::SU3vec(spinorfieldEvenOddElements, device));
			(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
		}
	}
protected:
	std::vector<const hardware::buffers::SU3vec *> spinorfields;
	const hardware::buffers::SU3vec * getOutSpinor() const
	{
		return spinorfields.back();
	}
};

struct EvenOddLinearCombinationTesterWithSquarenormAsResult : public EvenOddLinearCombinationTester
{
	EvenOddLinearCombinationTesterWithSquarenormAsResult(const std::string kernelName, const ParameterCollection parameterCollection, const EvenOddLinearCombinationTestParameters testParameters) :
		EvenOddLinearCombinationTester(kernelName, parameterCollection, testParameters) {};

	~EvenOddLinearCombinationTesterWithSquarenormAsResult()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinor());
	}
};

BOOST_AUTO_TEST_SUITE(SQUARENORM)

	struct SquarenormTestParameters: public NonEvenOddLinearCombinationTestParameters
	{
		SquarenormTestParameters(const LatticeExtents & latticeExtendsIn, const SpinorFillTypes & fillTypesIn) :
			NonEvenOddLinearCombinationTestParameters{calculateReferenceValues_squarenorm( getSpinorfieldSize(latticeExtendsIn), fillTypesIn),
			latticeExtendsIn, fillTypesIn} {};
	};

	class SquarenormTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
	  public:
		SquarenormTester(const ParameterCollection & parameterCollection, const SquarenormTestParameters & testParameters):
			NonEvenOddLinearCombinationTesterWithSquarenormAsResult("squarenorm", parameterCollection, testParameters) {
		}
	};

	BOOST_AUTO_TEST_CASE( SQUARENORM_1 )
	{
	   performTest<SquarenormTester, SquarenormTestParameters>( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_2 )
	{
	   performTest<SquarenormTester, SquarenormTestParameters>( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::ascendingComplex} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_1 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>( LatticeExtents{ns8, nt8}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_2 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>( LatticeExtents{ns12, nt12}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_3 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>( LatticeExtents{ns16, nt16}, SpinorFillTypes{ SpinorFillType::one} );
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	struct ScalarProductTestParameters: public NonEvenOddLinearCombinationTestParameters
	{
		ScalarProductTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn):
			NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_scalarProduct(getSpinorfieldSize(latticeExtentsIn),fillTypesIn), latticeExtentsIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2) {};
	};

	class ScalarProductTester: public NonEvenOddLinearCombinationTester{
	public:
	 ScalarProductTester(const ParameterCollection & parameterCollection, const ScalarProductTestParameters testParameters):
		 NonEvenOddLinearCombinationTester("scalar_product", parameterCollection, testParameters)
	{
			code->set_complex_to_scalar_product_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
			hmc_complex resultHost;
			complexNums.at(0)->dump(&resultHost);

			kernelResult[0] = resultHost.re;
			kernelResult[1] = resultHost.im;
	}
};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
		performTest<ScalarProductTester, ScalarProductTestParameters>( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		performTest<ScalarProductTester, ScalarProductTestParameters>( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex} );
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		performTest<ScalarProductTester, ScalarProductTestParameters>( LatticeExtents{ns8, nt8}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		performTest<ScalarProductTester, ScalarProductTestParameters>( LatticeExtents{ns12, nt12}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		performTest<ScalarProductTester, ScalarProductTestParameters>( LatticeExtents{ns16, nt16}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(ZERO)

struct ZeroTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	ZeroTestParameters(LatticeExtents latticeExtendsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_zero(), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

	class ZeroTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
	   public:
	ZeroTester(const ParameterCollection & parameterCollection, const ZeroTestParameters & testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("zero", parameterCollection, testParameters)
		{
			code->set_zero_spinorfield_device(getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		performTest<ZeroTester, ZeroTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(COLD)

struct ColdTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	ColdTestParameters(LatticeExtents latticeExtendsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_cold(false), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

	class ColdTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
	   public:
		ColdTester(const ParameterCollection & parameterCollection, const ColdTestParameters testParameters):
			NonEvenOddLinearCombinationTesterWithSquarenormAsResult("cold", parameterCollection, testParameters)
			{
				code->set_cold_spinorfield_device(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		performTest<ColdTester, ColdTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX)

	struct SaxTestParameters : public NonEvenOddLinearCombinationTestParameters
	{
		SaxTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_sax(getSpinorfieldSize(latticeExtentsIn), coefficientsIn), latticeExtentsIn, fillTypesIn, coefficientsIn, 2) {}
	};

	class SaxTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
	   public:
		SaxTester(const ParameterCollection & parameterCollection, const SaxTestParameters testParameters):
			NonEvenOddLinearCombinationTesterWithSquarenormAsResult("sax", parameterCollection, testParameters)
			{
				code->sax_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		performTest<SaxTester, SaxTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		performTest<SaxTester, SaxTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		performTest<SaxTester, SaxTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		performTest<SaxTester, SaxTestParameters> (LatticeExtents{ns16,nt4}, ComplexNumbers {{1.,1.}});
	}
	
	BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPY)

	struct SaxpyTestParameters: public NonEvenOddLinearCombinationTestParameters
	{
		SaxpyTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpy(getSpinorfieldSize(latticeExtentsIn), coefficientsIn), latticeExtentsIn, fillTypesIn, coefficientsIn, 3){}

	};

	class SaxpyTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
		public:
			SaxpyTester(const ParameterCollection & parameterCollection, const SaxpyTestParameters testParameters):
				NonEvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy", parameterCollection, testParameters)
				{
					code->saxpy_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
				}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (LatticeExtents{ns4, nt4}, ComplexNumbers{{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers{{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers{{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (LatticeExtents{ns8,nt8}, ComplexNumbers{{1.,1.}});

	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBYPZ)

	struct SaxpbypzTestParameters: public NonEvenOddLinearCombinationTestParameters
	{
		SaxpbypzTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpbypz(getSpinorfieldSize(latticeExtentsIn), coefficientsIn), latticeExtentsIn, fillTypesIn, coefficientsIn, 4){}
	};

	class SaxpbypzTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult{
		public:
			SaxpbypzTester(const ParameterCollection & parameterCollection, const SaxpbypzTestParameters testParameters):
				NonEvenOddLinearCombinationTesterWithSquarenormAsResult("saxpbypz", parameterCollection, testParameters)
		{
				code->saxpbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(3), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_1 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers{{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		performTest<SaxpbypzTester, SaxpbypzTestParameters> (LatticeExtents{ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	class GaussianTester: public SpinorStaggeredTester{
	   public:
		GaussianTester(std::string inputfile) : SpinorStaggeredTester("gaussian_spinorfield", inputfile, 1, 2){
			const hardware::buffers::Plain<su3vec> out(spinorfieldElements, device);
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);

			su3vec * outHost;
			outHost = new su3vec[spinorfieldElements * iterations];
			BOOST_REQUIRE(out);
				
			auto prng_buf = prng->get_buffers().at(0);
		
			double sum = 0;
			for (int i = 0; i < iterations; i++) {
				if(i%100==0) logger.info() << "Run kernel for the " << i << "th time";
				code->set_gaussian_spinorfield_device(&out, prng_buf);
				out.dump(&outHost[i * spinorfieldElements]);
				sum += count_sf(&outHost[i * spinorfieldElements], spinorfieldElements);
			}
			//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
			sum /= (iterations * spinorfieldElements * 6);
			kernelResult[0] = sum;

			if(calcVariance){
				double var = 0.;
				for (int i = 0; i < iterations; i++) {
				   var += calc_var_sf(&outHost[i * spinorfieldElements], spinorfieldElements, sum);
				}
				//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
				var /= (iterations * spinorfieldElements * 6);
				kernelResult[0] = sqrt(var);
			}
			
			/** 
			 * @TODO This piece of code contains actually tests for the RNG itself 
			 *       and should be moved elsewhere. 
			 */
			//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
			//(to big n and to small n are not good choices for this test)
			//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
			vector<vector<hmc_float>> samples;
			vector<hmc_float> tmp;
			vector<hmc_float> tmp2;
			for(int i=0; i<iterations; i++){
			  vector<hmc_float> tmp;
			  for(uint j=0; j<spinorfieldElements; j++){
			    tmp2=reals_from_su3vec(outHost[i*spinorfieldElements+j]);
			    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
			    tmp2.clear();
			  }
			  samples.push_back(tmp);
			  tmp.clear();
			}
			logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
			logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
			
			if(!calcVariance){
			  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
			}else{
			  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  variance_test_multiple_set(samples,2.,sqrt(0.5));
			  variance_test_multiple_set(samples,3.,sqrt(0.5));
			  variance_test_multiple_set(samples,4.,sqrt(0.5));
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
	    GaussianTester("gaussian_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_2 )
	{
	    GaussianTester("gaussian_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_3 )
	{
	    GaussianTester("gaussian_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_4 )
	{
	    GaussianTester("gaussian_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////           EVEN-ODD PRECONDITIONING             //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	struct ConvertEvenOddTestParameters: public SpinorStaggeredTestParameters
	{
		const bool fillEvenSites;
		ConvertEvenOddTestParameters(LatticeExtents latticeExtentsIn, const bool fillEvenSitesIn):
			SpinorStaggeredTestParameters(calculateReferenceValues_convert_eo(getEvenOddSpinorfieldSize(latticeExtentsIn), fillEvenSitesIn), latticeExtentsIn, true),
			fillEvenSites(fillEvenSitesIn) {}

	};

	class ConvertToEvenOddTester: public SpinorStaggeredTester{
		public:
			ConvertToEvenOddTester(const ParameterCollection & parameterCollection, ConvertEvenOddTestParameters testParameters):
				SpinorStaggeredTester("convert_to_eo", parameterCollection, testParameters)
		{
				const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
				const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);

				in.load(createSpinorfieldWithOnesAndZerosDependingOnSiteParity()); //createSpinorfieldWithOnesAndZerosDependingOnSiteParity should probably be changed as in SpinorTester
				code->convert_to_eoprec_device(&in2, &in3, &in);

				code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
				doubleBuffer->dump(&kernelResult[1]);
		}
};

	class ConvertEvenOddTester: public SpinorStaggeredTester{
	   public:
		ConvertEvenOddTester(std::string inputfile, bool from_eo) : SpinorStaggeredTester("convert fro/to eo", inputfile, 2){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);

			if(from_eo){
				in2.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity());
				in3.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity());
				code->convert_from_eoprec_device(&in2, &in3, &in);
				code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				referenceValue[0] = spinorfieldEvenOddElements*3;
				referenceValue[1] = 0.; //this has not actually a meaning
			}else{
				in.load(createSpinorfieldWithOnesAndZerosDependingOnSiteParity());
				code->convert_to_eoprec_device(&in2, &in3, &in);
				code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
				doubleBuffer->dump(&kernelResult[1]);
				
				if (evenOrOdd){
				  referenceValue[0] = spinorfieldEvenOddElements*3; 
				  referenceValue[1] = 0.;
				}else{ 
				  referenceValue[1] = spinorfieldEvenOddElements*3; 
				  referenceValue[0] = 0.;
				}
			}
		}
	};


	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
	    ConvertEvenOddTester("convert_eo_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
	    ConvertEvenOddTester("convert_eo_input_2", true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
	    ConvertEvenOddTester("convert_eo_input_1", false);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
	    ConvertEvenOddTester("convert_eo_input_2", false);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SQUARENORM_EO)

	struct SquarenormEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
	{
		SquarenormEvenOddTestParameters(const LatticeExtents & latticeExtendsIn, const SpinorFillTypes & fillTypesIn):
			EvenOddLinearCombinationTestParameters{calculateReferenceValues_squarenorm( getEvenOddSpinorfieldSize(latticeExtendsIn), fillTypesIn),latticeExtendsIn, fillTypesIn}{};
	};

	class SquarenormEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult{
		public:
			SquarenormEvenOddTester(const ParameterCollection & parameterCollection, const SquarenormEvenOddTestParameters & testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("squarenorm_eo", parameterCollection, testParameters){}
	};

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters>( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters>( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters>( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters>( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters>( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	struct ScalarProductEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
	{
		ScalarProductEvenOddTestParameters(const LatticeExtents & latticeExtentsIn, const SpinorFillTypes & fillTypesIn) :
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_scalarProduct( getEvenOddSpinorfieldSize(latticeExtentsIn), fillTypesIn), latticeExtentsIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2) {};
	};

	class ScalarProductEvenOddRealTester: public EvenOddLinearCombinationTester{
	public:
		ScalarProductEvenOddRealTester(const ParameterCollection & parameterCollection, const ScalarProductEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("scalar_product_real_eo", parameterCollection, testParameters)
		{
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);
			code->set_float_to_scalar_product_real_part_eoprec_device(spinorfields.at(0), spinorfields.at(1), &sqnorm);
			hmc_float resultTmp;
			sqnorm.dump(&resultTmp);
			kernelResult[0] = resultTmp;
			kernelResult[1] = 0;
		}
};

	class ScalarProductEvenOddComplexTester: public EvenOddLinearCombinationTester{
	public:
		ScalarProductEvenOddComplexTester(const ParameterCollection & parameterCollection, const ScalarProductEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("scalar_product_real_eo", parameterCollection, testParameters)
		{
			code->set_complex_to_scalar_product_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
			hmc_complex resultTmp;
			complexNums.at(0)->dump(&resultTmp);
			kernelResult[0] = resultTmp.re;
			kernelResult[1] = resultTmp.im;
		}
};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns8,nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns4,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_3 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_4 )
	{
		performTest<ScalarProductEvenOddComplexTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_1 )
	{
		performTest<ScalarProductEvenOddRealTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns8,nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_2 )
	{
		performTest<ScalarProductEvenOddRealTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_1 )
	{
		performTest<ScalarProductEvenOddRealTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_2 )
	{
		performTest<ScalarProductEvenOddRealTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_3 )
	{
		performTest<ScalarProductEvenOddRealTester, ScalarProductEvenOddTestParameters> ( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////
BOOST_AUTO_TEST_SUITE(ZERO_EO)

	struct ZeroEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
	{
		ZeroEvenOddTestParameters(const LatticeExtents latticeExtentsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_zero(), latticeExtentsIn, SpinorFillTypes{SpinorFillType::one}) {};
	};

	class ZeroEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		ZeroEvenOddTester(const ParameterCollection & parameterCollection, const ZeroEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("zero_eo", parameterCollection, testParameters)
			{
				code->set_zero_spinorfield_eoprec_device(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		performTest<ZeroEvenOddTester, ZeroEvenOddTestParameters> (LatticeExtents {ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_EO)

	struct ColdEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
	{
		ColdEvenOddTestParameters(const LatticeExtents latticeExtentsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_cold(true), latticeExtentsIn, SpinorFillTypes{SpinorFillType::one}) {};
	};

	class ColdEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult{
	   public:
		ColdEvenOddTester(const ParameterCollection & parameterCollection, const ColdEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("cold_eo", parameterCollection, testParameters)
		{
			code->set_cold_spinorfield_eoprec_device(getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		performTest<ColdEvenOddTester, ColdEvenOddTestParameters> (LatticeExtents{ns4,nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX_EO)

	struct SaxEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
	{
		SaxEvenOddTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_sax(getEvenOddSpinorfieldSize(latticeExtentsIn), coefficientsIn), latticeExtentsIn, fillTypesIn, coefficientsIn, 2){}
	};

	class SaxEvenOddComplexTester: EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
			SaxEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_eo_complex", parameterCollection, testParameters)
			{
				code->sax_eoprec_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
			}
	};

	class SaxArgEvenOddComplexTester: EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
			SaxArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_arg_eo_complex", parameterCollection, testParameters)
			{
				code->sax_eoprec_device(spinorfields.at(0), testParameters.coefficients.at(0), getOutSpinor());
			}
	};

	class SaxEvenOddRealTester: EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_eo_real", parameterCollection, testParameters)
		{
			hmc_complex complexNumbers;
			complexNums.at(0)->dump(&complexNumbers);
			code->sax_eoprec_device(spinorfields.at(0), complexNumbers.re, getOutSpinor());
		}
	};

	class SaxArgEvenOddRealTester: EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
			SaxArgEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_arg_eo_real", parameterCollection, testParameters)
		{
			code->sax_eoprec_device(spinorfields.at(0), testParameters.coefficients.at(0).re, getOutSpinor());
		}
	};

	class SaxVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
			SaxVecEvenOddTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_eo_vec", parameterCollection, testParameters)
		{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, testParameters.coefficients.at(0).re);
				const int index_alpha = 3;
				alpha_real_vec.load(&alpha_host_real_vec[0]);
				code->sax_eoprec_device(spinorfields.at(0), &alpha_real_vec, index_alpha, getOutSpinor());
		}

	};

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_1 )
	{
	    performTest<SaxEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_2 )
	{
	    performTest<SaxEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_3 )
	{
	    performTest<SaxEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_4 )
	{
	    performTest<SaxEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_1 )
	{
	    performTest<SaxArgEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_2 )
	{
	    performTest<SaxArgEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_3 )
	{
	    performTest<SaxArgEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_4 )
	{
	    performTest<SaxArgEvenOddComplexTester, SaxEvenOddTestParameters> (LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_1 )
	{
	    performTest<SaxEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_2 )
	{
	    performTest<SaxEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_3 )
	{
	    performTest<SaxEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_1 )
	{
	    performTest<SaxArgEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_2 )
	{
	    performTest<SaxArgEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_3 )
	{
	    performTest<SaxArgEvenOddRealTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_1 )
	{
	    performTest<SaxVecEvenOddTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers{{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_2 )
	{
	    performTest<SaxVecEvenOddTester, SaxEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_3 )
	{
	    performTest<SaxVecEvenOddTester, SaxEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	struct SaxpyEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
	{
		SaxpyEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpy(getEvenOddSpinorfieldSize(latticeExtentsIn), coefficientsIn),latticeExtentsIn, fillTypesIn, coefficientsIn, 3){}
	};

	class SaxpyEvenOddComplexTester: public EvenOddLinearCombinationTester
	{
	public:
		SaxpyEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("saxpy_eo_complex", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
		}
	};

	class SaxpyArgEvenOddComplexTester: public EvenOddLinearCombinationTester
	{
	public:
		SaxpyArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("saxpy_arg_eo_complex", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), getOutSpinor());
		}
	};

	class SaxpyEvenOddRealTester: public EvenOddLinearCombinationTester
	{
	public:
		SaxpyEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("saxpy_eo_real", parameterCollection, testParameters)
		{
			hmc_complex complexNumbers;
			complexNums.at(0)->dump(&complexNumbers);
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNumbers.re, getOutSpinor());
		}
	};

	class SaxpyArgEvenOddRealTester: public EvenOddLinearCombinationTester
	{
	public:
		SaxpyArgEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTester("saxpy_arg_eo_real", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0).re, getOutSpinor());
		}
	};

	class SaxpyVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
			SaxpyVecEvenOddTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
				EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_eo_vec", parameterCollection, testParameters)
		{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, testParameters.coefficients.at(0).re);
				const int index_alpha = 3;
				alpha_real_vec.load(&alpha_host_real_vec[0]);
				code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), &alpha_real_vec, index_alpha, getOutSpinor());
		}

	};

	class SaxpyEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpyEvenOddTester(std::string inputfile, int switcher):SpinorStaggeredTester("saxpy_eo or saxpy_arg_eo", inputfile, 1){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			
			if(switcher==0 || switcher==1){
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				alpha.load(&alpha_host);
				if(switcher==0)
				    code->saxpy_eoprec_device(&in, &in2, &alpha, &out);
				if(switcher==1)
				    code->saxpy_eoprec_device(&in, &in2, alpha_host, &out);
			}else if(switcher==2 || switcher==3){
				hardware::buffers::Plain<hmc_float> alpha_real(1, device);
				alpha_real.load(&alpha_host.re);
				if(switcher==2)
				    code->saxpy_eoprec_device(&in, &in2, &alpha_real, &out);
				if(switcher==3)
				    code->saxpy_eoprec_device(&in, &in2, alpha_host.re, &out);
			}else{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, alpha_host.re);
				const int index_alpha = 3;
				alpha_real_vec.load(&alpha_host_real_vec[0]); 
				code->saxpy_eoprec_device(&in,  &in2, &alpha_real_vec, index_alpha, &out);
			}
				
			calcSquarenormEvenOddAndStoreAsKernelResult(&out);

     /*
     print_staggeredfield_eo_to_textfile("ref_vec_saxpy1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
     logger.info() << "Produced the ref_vec_saxpy1_eo text file with the staggered field for the ref. code.";   
     print_staggeredfield_eo_to_textfile("ref_vec_saxpy2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
     logger.info() << "Produced the ref_vec_saxpy2_eo text file with the staggered field for the ref. code.";
     */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_3", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_4", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_5", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_7 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_7", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_8 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_8", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_9 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_9", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_10 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_11 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_12 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_12", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_13 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_13", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_14 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_14", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_15 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_16 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_16", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_17 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_17", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_18 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_18", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_3", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_4", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_5", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_7 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_7", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_8 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_8", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_9 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_9", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_10 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_11 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_12 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_12", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_13 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_13", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_14 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_14", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_15 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_16 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_16", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_17 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_17", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_18 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_18", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_3 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_6", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_4 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_10", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_5 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_11", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_6 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_15", 4);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBY_EO)

	class SaxpbyEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpbyEvenOddTester(std::string inputfile, int switcher):SpinorStaggeredTester("saxpy_eo or saxpy_arg_eo", inputfile, 1){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			
			if(switcher==0 || switcher==1){
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);
				alpha.load(&alpha_host);
				beta.load(&beta_host);
				if(switcher==0)
				    code->saxpby_eoprec_device(&in, &in2, &alpha, &beta, &out);
				if(switcher==1)
				    code->saxpby_eoprec_device(&in, &in2, alpha_host, beta_host, &out);
			}else if(switcher==2 || switcher==3){
				hardware::buffers::Plain<hmc_float> alpha_real(1, device);
				hardware::buffers::Plain<hmc_float> beta_real(1, device);
				alpha_real.load(&alpha_host.re);
				beta_real.load(&beta_host.re);
				if(switcher==2)
				    code->saxpby_eoprec_device(&in, &in2, &alpha_real, &beta_real, &out);
				if(switcher==3)
				    code->saxpby_eoprec_device(&in, &in2, alpha_host.re, beta_host.re, &out);
			}else{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, alpha_host.re);
				const int index_alpha = 3;
				hardware::buffers::Plain<hmc_float> beta_real_vec(5, device);
				std::vector<hmc_float> beta_host_real_vec(5, beta_host.re);
				const int index_beta = 2;
				alpha_real_vec.load(&alpha_host_real_vec[0]); 
				beta_real_vec.load(&beta_host_real_vec[0]); 
				code->saxpby_eoprec_device(&in, &in2, &alpha_real_vec, &beta_real_vec, index_alpha, index_beta, &out);
			}
			
			calcSquarenormEvenOddAndStoreAsKernelResult(&out);

     /*
     print_staggeredfield_eo_to_textfile("ref_vec_saxpby1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
     logger.info() << "Produced the ref_vec_saxpby1_eo text file with the staggered field for the ref. code.";   
     print_staggeredfield_eo_to_textfile("ref_vec_saxpby2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
     logger.info() << "Produced the ref_vec_saxpby2_eo text file with the staggered field for the ref. code.";
     */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBY_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_2", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_4", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_5", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_6", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_7", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_8", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_9 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_10 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_10", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_11 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_12 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_12", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_13 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_13", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_14 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_14", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_15 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_15", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_16 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_16", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_17 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_17", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_18 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_19 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_19", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_20 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_21 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_21", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_22 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_22", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_23 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_23", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_24 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_24", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_25 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_25", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_26 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_27 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_27", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_28 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_29 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_29", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_30 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_30", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_31 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_31", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_32 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_32", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_33 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_33", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_34 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_34", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_2", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_4", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_5", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_6", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_7", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_8", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_9 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_10 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_10", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_11 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_12 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_12", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_13 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_13", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_14 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_14", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_15 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_15", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_16 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_16", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_17 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_17", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_18 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_19 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_19", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_20 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_21 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_21", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_22 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_22", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_23 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_23", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_24 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_24", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_25 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_25", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_26 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_27 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_27", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_28 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_29 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_29", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_30 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_30", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_31 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_31", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_32 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_32", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_33 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_33", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_34 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_34", 1);
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 4);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBYPZ_EO)

	class SaxpbypzEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpbypzEvenOddTester(std::string inputfile, bool switcher=true) : 
		                               SpinorStaggeredTester("saxsbypz_eo", inputfile){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);
			hardware::buffers::Plain<hmc_complex> beta(1, device);

			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			in3.load(createSpinorfield(spinorfieldEvenOddElements, 789));
			alpha.load(&alpha_host);
			beta.load(&beta_host);

			(switcher) ? code->saxpbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out) :
			             code->saxpbypz_eoprec_device(&in, &in2, &in3, alpha_host, beta_host, &out);

			calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			
  /*
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
  logger.info() << "Produced the ref_vec_saxpbypz1_eo text file with the staggered field for the ref. code."; 
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
  logger.info() << "Produced the ref_vec_saxpbypz2_eo text file with the staggered field for the ref. code.";  
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz3_eo", createSpinorfield(spinorfieldEvenOddElements, 789)); 
  logger.info() << "Produced the ref_vec_saxpbypz3_eo text file with the staggered field for the ref. code.";
  */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_1 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_2 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_3 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_4 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_5 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_6 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_7 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_8 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_9 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_10 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_10");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_11 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_11");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_12 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_12");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_13 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_13");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_14 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_14");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_15 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_15");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_16 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_16");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_17 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_17");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_18 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_18");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_19 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_19");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_20 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_20");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_21 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_21");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_22 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_22");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_23 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_23");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_24 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_24");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_25 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_25");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_26 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_26");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_27 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_27");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_28 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_28");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_29 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_29");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_30 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_30");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_31 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_31");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_32 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_32");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_33 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_33");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_34 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_34");
	}

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_1 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_1", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_2 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_2", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_3 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_3", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_4 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_4", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_5 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_5", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_6 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_6", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_7 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_7", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_8 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_8", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_9 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_9", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_10 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_10", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_11 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_11", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_12 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_12", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_13 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_13", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_14 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_14", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_15 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_15", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_16 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_16", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_17 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_17", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_18 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_18", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_19 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_19", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_20 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_20", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_21 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_21", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_22 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_22", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_23 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_23", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_24 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_24", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_25 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_25", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_26 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_26", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_27 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_27", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_28 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_28", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_29 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_29", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_30 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_30", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_31 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_31", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_32 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_32", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_33 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_33", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_34 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_34", false);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	class GaussianEvenOddTester: public SpinorStaggeredTester{
	   public:
		GaussianEvenOddTester(std::string inputfile) : SpinorStaggeredTester("gaussian_spinorfield_eo", inputfile, 1, 2){
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);

			su3vec * outHost;
			outHost = new su3vec[spinorfieldEvenOddElements * iterations];
			BOOST_REQUIRE(out);
				
			auto prng_buf = prng->get_buffers().at(0);
		
			double sum = 0;
			for (int i = 0; i < iterations; i++) {
				if(i%200==0) logger.info() << "Run kernel for the " << i << "th time";
				code->set_gaussian_spinorfield_eoprec_device(&out, prng_buf);
				out.dump(&outHost[i * spinorfieldEvenOddElements]);
				sum += count_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements);
			}
			//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
			sum /= (iterations * spinorfieldEvenOddElements * 6);
			kernelResult[0] = sum;

			if(calcVariance){
				double var = 0.;
				for (int i = 0; i < iterations; i++) {
				   var += calc_var_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements, sum);
				}
				//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
				var /= (iterations * spinorfieldEvenOddElements * 6);
				kernelResult[0] = sqrt(var);
			}
			
			/** 
			 * @TODO This piece of code contains actually tests for the RNG itself 
			 *       and should be moved elsewhere. 
			 */
			//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
			//(to big n and to small n are not good choices for this test)
			//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
			vector<vector<hmc_float>> samples;
			vector<hmc_float> tmp;
			vector<hmc_float> tmp2;
			for(int i=0; i<iterations; i++){
			  vector<hmc_float> tmp;
			  for(uint j=0; j<spinorfieldEvenOddElements; j++){
			    tmp2=reals_from_su3vec(outHost[i*spinorfieldEvenOddElements+j]);
			    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
			    tmp2.clear();
			  }
			  samples.push_back(tmp);
			  tmp.clear();
			}
			logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
			logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
			
			if(!calcVariance){
			  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
			}else{
			  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  variance_test_multiple_set(samples,2.,sqrt(0.5));
			  variance_test_multiple_set(samples,3.,sqrt(0.5));
			  variance_test_multiple_set(samples,4.,sqrt(0.5));
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_2 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_3 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_4 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX_VEC_AND_SQNORM)

	class SaxVecAndSqnormEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxVecAndSqnormEvenOddTester(std::string inputfile) : SpinorStaggeredTester("sax_vectorized_and_squarenorm_eoprec",inputfile, 1){
			int NUM_EQS = parameters->get_md_approx_ord();
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements*NUM_EQS, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements*NUM_EQS, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			
			hardware::buffers::Plain<hmc_float> sqnorm(NUM_EQS, device);
			hardware::buffers::Plain<hmc_float> alpha(NUM_EQS, device);
			std::vector<hmc_float> alpha_vec_host(NUM_EQS, alpha_host.re);
			for(uint i=0; i<alpha_vec_host.size(); i++)
			    alpha_vec_host[i] += i*alpha_host.im;
			alpha.load(&alpha_vec_host[0]);
			
			code->sax_vectorized_and_squarenorm_eoprec_device(&in, &alpha, NUM_EQS, &sqnorm);
			std::vector<hmc_float> cpu_res(NUM_EQS);
			sqnorm.dump(&cpu_res[0]);

			kernelResult[0] = std::accumulate(cpu_res.begin(), cpu_res.end(), 0.0);

  /*
  print_staggeredfield_eo_to_textfile("ref_vec_sax_and_sq_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
  logger.info() << "Produced the ref_vec_sax_and_sq_eo text file with the staggered field for the ref. code.";
  */
		}
	};

	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_1 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_2 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_3 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_4 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_5 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_6 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_7 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_8 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_9 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_10 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_10");
	}

BOOST_AUTO_TEST_SUITE_END()
