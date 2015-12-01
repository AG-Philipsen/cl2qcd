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

const ReferenceValues calculateReferenceValue_saxpby(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{alphaIn.at(0).re + alphaIn.at(1).re, alphaIn.at(0).im + alphaIn.at(1).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_saxpbypz(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{1. + alphaIn.at(0).re + alphaIn.at(1).re, 0. + alphaIn.at(0).im + alphaIn.at(1).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_gaussian()
{
	return ReferenceValues{1e-3, .71};
}

const ReferenceValues calculateReferenceValues_convert_eo(const int latticeVolume, const bool fillEvenSites)
{
	double nonTrivialValue = latticeVolume * 3.;
	return ReferenceValues {fillEvenSites ? nonTrivialValue : 0, fillEvenSites ? 0. : nonTrivialValue};
}

const ReferenceValues calculateReferenceValues_convertFromEvenOdd(const int latticeVolume)
{
	return ReferenceValues {latticeVolume * 3. / 2.};
}

const ReferenceValues calculateReferenceValue_sax_vec_and_sqnorm(const int latticeVolume, const ComplexNumbers alphaIn, int numEqsIn)
{
	std::vector<double> alpha;
	for (int i=0; i<numEqsIn; i++)	alpha.push_back((alphaIn.at(0).re + i * alphaIn.at(0).im)*(alphaIn.at(0).re + i * alphaIn.at(0).im));
	return ReferenceValues {latticeVolume * 3. * std::accumulate(alpha.begin(), alpha.end(), 0.0)};
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

template<typename TesterClass, typename ParameterClass> void performTest(LatticeExtents latticeExtendsIn, const ComplexNumbers alphaIn, const int numEqsIn )
{
	ParameterClass parametersForThisTest(latticeExtendsIn, alphaIn, numEqsIn);
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
		const hardware::HardwareParametersMockup hardwareParameters(4,4, true);//this implicitly sets useEvenOdd(false) if true is not there
		const hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(4,4);
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

struct GaussianTestParameters: public NonEvenOddLinearCombinationTestParameters
{
	GaussianTestParameters(const LatticeExtents latticeExtentsIn, const ComparisonType & typeOfComparisonIn) :
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_gaussian(), latticeExtentsIn, 1, typeOfComparisonIn),
		iterations(4000){};

	const unsigned int iterations;
};

struct PrngTester: public NonEvenOddLinearCombinationTester
{
	PrngTester(const std::string kernelName, const ParameterCollection parameterCollection, const GaussianTestParameters & testParameters):
		NonEvenOddLinearCombinationTester(kernelName, parameterCollection, testParameters),
				hostSeed( parameterCollection.kernelParameters.getHostSeed() ),
				useSameRandomNumbers(parameterCollection.hardwareParameters.useSameRandomNumbers())
	{
		prng_init(hostSeed);
		prngStates = new hardware::buffers::PRNGBuffer(device, useSameRandomNumbers );
		auto codePrng = device->getPrngCode();
		codePrng->initialize(prngStates, hostSeed);
	}
	~PrngTester()
	{
		delete prngStates;
	}
protected:
	const hardware::buffers::PRNGBuffer* prngStates;
private:
	uint32_t hostSeed;
	bool useSameRandomNumbers;
};

class GaussianTester: public PrngTester{
public:
	GaussianTester(std::string kernelName, const ParameterCollection & parameterCollection, const GaussianTestParameters & testParameters, const int numberOfElements):
		PrngTester(kernelName, parameterCollection, testParameters),
		numberOfElements(numberOfElements), mean(0.), variance(0.),
		hostOutput(std::vector<su3vec>(numberOfElements * testParameters.iterations)), testParameters(testParameters){}

	~GaussianTester()
	{
		calculateMean();
		calculateVariance();
		KolmogorovSmirnov();

		kernelResult[0] = mean;
		kernelResult[1] = sqrt(variance);
	}
	double normalize(double valueIn)
	{
		return valueIn/= testParameters.iterations * numberOfElements * 6;
	}

	void calculateMean()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
//				if(i%100==0) logger.info() << "Run kernel for the " << i << "th time";
			mean += count_sf(&hostOutput[i * numberOfElements], numberOfElements);
		}
		mean = normalize(mean);
	}

	void calculateVariance()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			variance += calc_var_sf(&hostOutput[i * numberOfElements], numberOfElements, mean);
		}
		variance = normalize(variance);
	}

	void KolmogorovSmirnov()
	{
		vector<vector<hmc_float>> samples;
		vector<hmc_float> tmp;
		vector<hmc_float> tmp2;
		for(unsigned int i=0; i<testParameters.iterations; i++){
		  vector<hmc_float> tmp;
		  for(int j=0; j<numberOfElements; j++){
		    tmp2=reals_from_su3vec(hostOutput[i*numberOfElements+j]);
		    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
		    tmp2.clear();
		  }
		  samples.push_back(tmp);
		  tmp.clear();
		}
		logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
		logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";

		//Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
		//Note that in the test BOOST_CHECK is used.
		mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
		mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
		mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
		//Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
		//Note that in the test BOOST_CHECK is used.
		variance_test_multiple_set(samples,2.,sqrt(0.5));
		variance_test_multiple_set(samples,3.,sqrt(0.5));
		variance_test_multiple_set(samples,4.,sqrt(0.5));
	}
	protected:
		const int numberOfElements;
		double mean, variance;
		std::vector<su3vec> hostOutput;
		const GaussianTestParameters & testParameters;
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

	struct NonEvenOddGaussianStaggeredSpinorfieldTester: public GaussianTester
	{
		NonEvenOddGaussianStaggeredSpinorfieldTester(const ParameterCollection & parameterCollection, const GaussianTestParameters testParameters):
			GaussianTester("generate_gaussian_staggeredspinorfield", parameterCollection, testParameters, testParameters.getSpinorfieldSize()) {}
		~NonEvenOddGaussianStaggeredSpinorfieldTester()
		{
			const hardware::buffers::Plain<su3vec> outSpinor(spinorfieldElements, device);
			for (unsigned int i = 0; i < testParameters.iterations; i++){
				code->set_gaussian_spinorfield_device(&outSpinor,prngStates);
				outSpinor.dump(&hostOutput[i * numberOfElements]);
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
	    performTest<NonEvenOddGaussianStaggeredSpinorfieldTester, GaussianTestParameters>(LatticeExtents{ns4,nt4}, ComparisonType::smallerThan);
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

	class ConvertToEvenOddTester: public SpinorStaggeredTester
	{
		public:
			ConvertToEvenOddTester(const ParameterCollection & parameterCollection, ConvertEvenOddTestParameters testParameters):
				SpinorStaggeredTester("convert_to_eo", parameterCollection, testParameters)
		{
				const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
				const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);

				in.load(createSpinorfieldWithOnesAndZerosDependingOnSiteParity( testParameters.fillEvenSites ));
				code->convert_to_eoprec_device(&in2, &in3, &in);

				code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
				doubleBuffer->dump(&kernelResult[1]);
		}
	};

	struct ConvertFromEvenOddTestParameters: public SpinorStaggeredTestParameters
	{
		const bool fillEvenSites;
		ConvertFromEvenOddTestParameters(LatticeExtents latticeExtentsIn, const bool fillEvenSitesIn):
			SpinorStaggeredTestParameters(calculateReferenceValues_convertFromEvenOdd(getSpinorfieldSize(latticeExtentsIn)), latticeExtentsIn, true),
			fillEvenSites(fillEvenSitesIn) {}

	};

	class ConvertFromEvenOddTester: public SpinorStaggeredTester
	{
	public:
		ConvertFromEvenOddTester(const ParameterCollection & parameterCollection, ConvertFromEvenOddTestParameters testParameters):
			SpinorStaggeredTester("convert_from_eo", parameterCollection, testParameters)
		{
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);

			in2.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity( testParameters.fillEvenSites ));
			in3.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity( testParameters.fillEvenSites ));
			code->convert_from_eoprec_device(&in2, &in3, &in);
			code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);

		}
	};

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
	    performTest<ConvertToEvenOddTester, ConvertEvenOddTestParameters> (LatticeExtents {ns4, nt4}, true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
	    performTest<ConvertToEvenOddTester, ConvertEvenOddTestParameters> (LatticeExtents {ns8, nt4}, false);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		performTest<ConvertFromEvenOddTester, ConvertFromEvenOddTestParameters>(LatticeExtents{ns4, nt4}, true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
		performTest<ConvertFromEvenOddTester, ConvertFromEvenOddTestParameters>(LatticeExtents{ns8, nt4}, false);
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
		SaxpyEvenOddTestParameters(LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpy(getEvenOddSpinorfieldSize(latticeExtentsIn), coefficientsIn),latticeExtentsIn, fillTypesIn, coefficientsIn, 3){}
	};

	class SaxpyEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpyEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_eo_complex", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
		}
	};

	class SaxpyArgEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpyArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_arg_eo_complex", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), getOutSpinor());
		}
	};

	class SaxpyEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpyEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_eo_real", parameterCollection, testParameters)
		{
			hmc_complex complexNumbers;
			complexNums.at(0)->dump(&complexNumbers);
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNumbers.re, getOutSpinor());
		}
	};

	class SaxpyArgEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpyArgEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_arg_eo_real", parameterCollection, testParameters)
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

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_1 )
	{
	    performTest<SaxpyEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_2 )
	{
	    performTest<SaxpyEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_3 )
	{
	    performTest<SaxpyEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_4 )
	{
	    performTest<SaxpyEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_1 )
	{
	    performTest<SaxpyArgEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_2 )
	{
	    performTest<SaxpyArgEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_3 )
	{
	    performTest<SaxpyArgEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_4 )
	{
	    performTest<SaxpyArgEvenOddComplexTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_1 )
	{
	    performTest<SaxpyEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_2 )
	{
	    performTest<SaxpyEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_3 )
	{
	    performTest<SaxpyEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_1 )
	{
	    performTest<SaxpyArgEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_2 )
	{
	    performTest<SaxpyArgEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_3 )
	{
	    performTest<SaxpyArgEvenOddRealTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_1 )
	{
	    performTest<SaxpyVecEvenOddTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_2 )
	{
	    performTest<SaxpyVecEvenOddTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_3 )
	{
	    performTest<SaxpyVecEvenOddTester, SaxpyEvenOddTestParameters> (LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBY_EO)

	struct SaxpbyEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
	{
		SaxpbyEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValue_saxpby(getEvenOddSpinorfieldSize(latticeExtentsIn), coefficientsIn),latticeExtentsIn, fillTypesIn, coefficientsIn, 3 ) {}
	};

	class SaxpbyEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbyEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpbyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpby_eo_complex", parameterCollection, testParameters)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
	};

	class SaxpbyArgEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbyArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const SaxpbyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpby_arg_eo_complex", parameterCollection, testParameters)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), testParameters.coefficients.at(1), getOutSpinor());
		}
	};

	class SaxpbyEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbyEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpbyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpby_eo_real", parameterCollection, testParameters)
		{
			hmc_complex complexNumbers1;
			hmc_complex complexNumbers2;
			complexNums.at(0)->dump(&complexNumbers1);
			complexNums.at(1)->dump(&complexNumbers2);
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNumbers1.re, complexNumbers2.re, getOutSpinor());
		}
	};

	class SaxpbyArgEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbyArgEvenOddRealTester(const ParameterCollection & parameterCollection, const SaxpbyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpby_arg_eo_real", parameterCollection, testParameters)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0).re, testParameters.coefficients.at(1).re, getOutSpinor());
		}
	};

	class SaxpbyVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbyVecEvenOddTester(const ParameterCollection parameterCollection, const SaxpbyEvenOddTestParameters testParameters):
			EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpby_eo_vec", parameterCollection, testParameters)
		{
			hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
			std::vector<hmc_float> alpha_host_real_vec(5, testParameters.coefficients.at(0).re);
			const int index_alpha = 3;
			hardware::buffers::Plain<hmc_float> beta_real_vec(5, device);
			std::vector<hmc_float> beta_host_real_vec(5, testParameters.coefficients.at(1).re);
			const int index_beta = 2;
			alpha_real_vec.load(&alpha_host_real_vec[0]);
			beta_real_vec.load(&beta_host_real_vec[0]);
			code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), &alpha_real_vec, &beta_real_vec, index_alpha, index_beta, getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_1 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_2 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_3 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_4 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_5 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_6 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_7 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_8 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_9 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_10 )
	{
	    performTest<SaxpbyEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_1 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_2 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_3 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_4 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_5 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_6 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_7 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_8 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_9 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_10 )
	{
	    performTest<SaxpbyArgEvenOddComplexTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_1 )
	{
	    performTest<SaxpbyEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_2 )
	{
	    performTest<SaxpbyEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_3 )
	{
	    performTest<SaxpbyEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_4 )
	{
	    performTest<SaxpbyEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_5 )
	{
	    performTest<SaxpbyEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_1 )
	{
	    performTest<SaxpbyArgEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_2 )
	{
	    performTest<SaxpbyArgEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_3 )
	{
	    performTest<SaxpbyArgEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_4 )
	{
	    performTest<SaxpbyArgEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_5 )
	{
	    performTest<SaxpbyArgEvenOddRealTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_1 )
	{
	    performTest<SaxpbyVecEvenOddTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_2 )
	{
	    performTest<SaxpbyVecEvenOddTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_3 )
	{
	    performTest<SaxpbyVecEvenOddTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_4 )
	{
	    performTest<SaxpbyVecEvenOddTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_5 )
	{
	    performTest<SaxpbyVecEvenOddTester,SaxpbyEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBYPZ_EO)

	struct SaxpbypzEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
	{
		SaxpbypzEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			EvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpbypz(getEvenOddSpinorfieldSize(latticeExtentsIn),coefficientsIn), latticeExtentsIn, fillTypesIn, coefficientsIn, 4) {}
	};

	class SaxpbypzEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbypzEvenOddTester(const ParameterCollection & parameterCollection, SaxpbypzEvenOddTestParameters testParameters):
			 EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpbypz_eo", parameterCollection, testParameters) {
				code->saxpbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}

	};

	class SaxpbypzArgEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
	{
	public:
		SaxpbypzArgEvenOddTester(const ParameterCollection & parameterCollection, SaxpbypzEvenOddTestParameters testParameters):
			 EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpbypz_arg_eo", parameterCollection, testParameters) {
				code->saxpbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), testParameters.coefficients.at(0), testParameters.coefficients.at(1), getOutSpinor());
		}

	};

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_1 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_2 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_3 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_4 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_5 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_6 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_7 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_8 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_9 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_10 )
	{
	    performTest<SaxpbypzEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_1 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_2 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_3 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_4 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_5 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_6 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_7 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_8 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_9 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_10 )
	{
	    performTest<SaxpbypzArgEvenOddTester,SaxpbypzEvenOddTestParameters> (LatticeExtents {ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	struct EvenOddGaussianStaggeredSpinorfieldTester: public GaussianTester
	{
		EvenOddGaussianStaggeredSpinorfieldTester(const ParameterCollection & parameterCollection, const GaussianTestParameters testParameters):
			GaussianTester("generate_gaussian_staggeredspinorfield", parameterCollection, testParameters, testParameters.getSpinorfieldSize()) {}
		~EvenOddGaussianStaggeredSpinorfieldTester()
		{
			const hardware::buffers::Plain<su3vec> outSpinor(spinorfieldElements, device);
			for (unsigned int i = 0; i < testParameters.iterations; i++){
				code->set_gaussian_spinorfield_device(&outSpinor,prngStates);
				outSpinor.dump(&hostOutput[i * numberOfElements]);
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
	    performTest<EvenOddGaussianStaggeredSpinorfieldTester, GaussianTestParameters>(LatticeExtents{ns4,nt4}, ComparisonType::smallerThan);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX_VEC_AND_SQNORM)

	struct SaxVecAndSqnormEvenOddTestParameters: public EvenOddSpinorStaggeredTestParameters
	{
		SaxVecAndSqnormEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const ComplexNumbers coefficientsIn, const int numEqsIn):
		    EvenOddSpinorStaggeredTestParameters(calculateReferenceValue_sax_vec_and_sqnorm(getEvenOddSpinorfieldSize(latticeExtentsIn), coefficientsIn, numEqsIn), latticeExtentsIn),
			coefficients(coefficientsIn),
			numEqs(numEqsIn){};
			ComplexNumbers coefficients;
			const int numEqs;
	};

	class SaxVecAndSqnormEvenOddTester: public SpinorStaggeredTester
	{
	public:
		SaxVecAndSqnormEvenOddTester(const ParameterCollection & parameterCollection, const SaxVecAndSqnormEvenOddTestParameters testParameters):
			SpinorStaggeredTester("sax_vectorized_and_squarenorm_eoprec", parameterCollection, testParameters)
		{
			int NUM_EQS = testParameters.numEqs;
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements*NUM_EQS, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements*NUM_EQS, device);
			in.load(createSpinorfield(SpinorFillType::one));

			hardware::buffers::Plain<hmc_float> sqnorm(NUM_EQS, device);
			hardware::buffers::Plain<hmc_float> alpha(NUM_EQS, device);
			std::vector<hmc_float> alpha_vec_host(NUM_EQS, testParameters.coefficients.at(0).re);
			for(uint i=0; i<alpha_vec_host.size(); i++)
			    alpha_vec_host[i] += i*testParameters.coefficients.at(0).im;
			alpha.load(&alpha_vec_host[0]);

			code->sax_vectorized_and_squarenorm_eoprec_device(&in, &alpha, NUM_EQS, &sqnorm);
			std::vector<hmc_float> cpu_res(NUM_EQS);
			sqnorm.dump(&cpu_res[0]);

			kernelResult[0] = std::accumulate(cpu_res.begin(), cpu_res.end(), 0.0);

		}
	};

	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_1 )
	{
		performTest<SaxVecAndSqnormEvenOddTester, SaxVecAndSqnormEvenOddTestParameters> (LatticeExtents {ns4,nt4}, ComplexNumbers {{0.,0.}}, 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_2 )
	{
		performTest<SaxVecAndSqnormEvenOddTester, SaxVecAndSqnormEvenOddTestParameters> (LatticeExtents {ns4,nt4}, ComplexNumbers {{1.,0.}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_3 )
	{
		performTest<SaxVecAndSqnormEvenOddTester, SaxVecAndSqnormEvenOddTestParameters> (LatticeExtents {ns4,nt4}, ComplexNumbers {{1.,0.1}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_4 )
	{
		performTest<SaxVecAndSqnormEvenOddTester, SaxVecAndSqnormEvenOddTestParameters> (LatticeExtents {ns4,nt4}, ComplexNumbers {{0.,1.}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_5 )
	{
		performTest<SaxVecAndSqnormEvenOddTester, SaxVecAndSqnormEvenOddTestParameters> (LatticeExtents {ns4,nt4}, ComplexNumbers {{0.01,0.025}}, 15);
	}

BOOST_AUTO_TEST_SUITE_END()
