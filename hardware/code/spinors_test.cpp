/*
 * Copyright 2012, 2013, 2014, 2015 Christopher Pinke, Matthias Bach,
 *         Francesca Cuteri
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
#define BOOST_TEST_MODULE HARDWARE_CODE_SPINORS

#include "SpinorTester.hpp"
#include "../../host_functionality/host_random.h" //@todo: remove this in the end!

const ReferenceValues calculateReferenceValues_globalSquarenorm(const int latticeVolume, const SpinorFillTypes fillTypesIn)
{
	switch( fillTypesIn.at(0) )
	{
		case SpinorFillType::one :
		{
			return ReferenceValues{latticeVolume * 12.};
		}
		case SpinorFillType::ascendingComplex:
		{
			return ReferenceValues{latticeVolume * sumOfIntegersSquared(24)};
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
		return ReferenceValues{latticeVolume * 12., 0.};
	}
	else if( fillTypesIn.at(0) == SpinorFillType::one and fillTypesIn.at(1) == SpinorFillType::ascendingComplex )
	{
		return ReferenceValues{latticeVolume * sumOfIntegers(1,23,2), latticeVolume * sumOfIntegers(2,24,2)};
	}
	else if( fillTypesIn.at(0) == SpinorFillType::ascendingComplex and fillTypesIn.at(1) == SpinorFillType::one )
	{
		return ReferenceValues{latticeVolume * sumOfIntegers(1,23,2), -latticeVolume * sumOfIntegers(2,24,2)};
	}
	else if ( fillTypesIn.at(0) == SpinorFillType::ascendingComplex and fillTypesIn.at(1) == SpinorFillType::ascendingComplex  )
	{
		return ReferenceValues{latticeVolume * sumOfIntegersSquared(24), 0.};
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
	return ReferenceValues{ (alphaIn.at(0).im * alphaIn.at(0).im + alphaIn.at(0).re * alphaIn.at(0).re) * latticeVolume * sumOfIntegersSquared(24)};
}

const ReferenceValues calculateReferenceValues_saxpy(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{1. - alphaIn.at(0).re, 0. - alphaIn.at(0).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_saxsbypz(const int latticeVolume, const ComplexNumbers alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, ComplexNumbers {{1. + alphaIn.at(0).re + alphaIn.at(1).re, 0. + alphaIn.at(0).im + alphaIn.at(1).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_gaussian()
{
	return ReferenceValues{ 1e-3, .5};
}

const ReferenceValues calculateReferenceValues_convert_eo(const int latticeVolume, const bool fillEvenSites)
{
	double nonTrivialValue = latticeVolume * 12.;
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

struct NonEvenOddLinearCombinationTestParameters : public NonEvenOddSpinorTestParameters, LinearCombinationTestParameters
{
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn) :
		NonEvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn), LinearCombinationTestParameters{coefficientsIn, numberOfSpinorsIn} {};
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const size_t numberOfSpinorsIn, const ComparisonType typeOfComparisonIn):
		NonEvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisonIn), LinearCombinationTestParameters{numberOfSpinorsIn}{};
	NonEvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		NonEvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn) {};
};

struct EvenOddLinearCombinationTestParameters : public EvenOddSpinorTestParameters, LinearCombinationTestParameters
{
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn) :
		EvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn), LinearCombinationTestParameters{coefficientsIn, numberOfSpinorsIn} {};
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const size_t numberOfSpinorsIn, const ComparisonType typeOfComparisonIn):
		EvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisonIn), LinearCombinationTestParameters{numberOfSpinorsIn}{};
	EvenOddLinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		EvenOddSpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn) {};
};

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.needEvenOdd);
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

//@todo: this should take the coefficients from another class or be a template class!
struct NonEvenOddLinearCombinationTester3 : public NonEvenOddSpinorTester
{
	NonEvenOddLinearCombinationTester3(const std::string kernelName, const ParameterCollection pC, const NonEvenOddLinearCombinationTestParameters tP):
		NonEvenOddSpinorTester(kernelName, pC, tP)
	{
		loadCoefficients(tP);
		for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::Plain<spinor>(elements, device));
			(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(tP.fillTypes.at(number)));
		}
	}
	~NonEvenOddLinearCombinationTester3()
	{
		calcSquarenormAndStoreAsKernelResult(getOutSpinor());
	}
	protected:
		std::vector<const hardware::buffers::Plain<spinor> *> spinorfields;
		std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
		const hardware::buffers::Plain<spinor> * getOutSpinor() const
		{
			return spinorfields.back();
		}
	private:
		void loadCoefficients(const LinearCombinationTestParameters & testParameters)
		{
			for (auto coefficient : testParameters.coefficients)
			{
				complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, SpinorTester::device));
				complexNums.back()->load(&coefficient);
			}
		}
};

struct EvenOddLinearCombinationTester3 : public EvenOddSpinorTester
{
	EvenOddLinearCombinationTester3(const std::string kernelName, const ParameterCollection pC, const EvenOddLinearCombinationTestParameters tP):
		EvenOddSpinorTester(kernelName, pC, tP)
	{
		loadCoefficients(tP);
		for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::Spinor(elements, device));
			(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(tP.fillTypes.at(number)));
		}
	}
	~EvenOddLinearCombinationTester3()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinor());
	}
	protected:
		std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
		std::vector<const hardware::buffers::Spinor *> spinorfields;
		const hardware::buffers::Spinor * getOutSpinor() const
		{
			return spinorfields.back();
		}
	private:
		void loadCoefficients(const LinearCombinationTestParameters & testParameters)
		{
			for (auto coefficient : testParameters.coefficients)
			{
				complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, SpinorTester::device));
				complexNums.back()->load(&coefficient);
			}
		}
};

struct LinearCombinationTester: public SpinorTester
{
	LinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const NonEvenOddLinearCombinationTestParameters testParameters):
		SpinorTester(kernelName, parameterCollection, testParameters)
	{
		loadCoefficients(testParameters);
	}
	LinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const EvenOddLinearCombinationTestParameters testParameters):
		SpinorTester(kernelName, parameterCollection, testParameters)
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

struct NonEvenOddLinearCombinationTester: public LinearCombinationTester
{
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection parameterCollection, const NonEvenOddLinearCombinationTestParameters testParameters):
		LinearCombinationTester(kernelName, parameterCollection, testParameters)
	{
		for( size_t number = 0; number < testParameters.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::Plain<spinor>(spinorfieldElements, device));
			(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
		}
	}
protected:
	std::vector<const hardware::buffers::Plain<spinor> *> spinorfields;
	const hardware::buffers::Plain<spinor> * getOutSpinor() const
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

struct EvenOddLinearCombinationTester: public LinearCombinationTester
{
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection & parameterCollection, const EvenOddLinearCombinationTestParameters & testParameters):
		LinearCombinationTester(kernelName, parameterCollection, testParameters)
	{
		for( size_t number = 0; number < testParameters.numberOfSpinors ; number ++)
		{
			spinorfields.push_back(new hardware::buffers::Spinor(spinorfieldEvenOddElements, device));
			(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
		}
	}
protected:
	std::vector<const hardware::buffers::Spinor *> spinorfields;
	const hardware::buffers::Spinor * getOutSpinor() const
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

struct PrngTester: public SpinorTester
{
	PrngTester(const std::string kernelName, const ParameterCollection parameterCollection, const SpinorTestParameters & testParameters):
				SpinorTester(kernelName, parameterCollection, testParameters),
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


struct GaussianTestParameters : public SpinorTestParameters
{
	GaussianTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const ComparisonType & typeOfComparisonIn, const bool needsEvenOdd) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisonIn, needsEvenOdd), iterations(100) {};

	const unsigned int iterations;
};

struct GaussianTester: public PrngTester
{
	GaussianTester(const std::string kernelName, const ParameterCollection parameterCollection, const GaussianTestParameters & testParameters, const int numberOfElements):
				PrngTester(kernelName, parameterCollection, testParameters),
				numberOfElements(numberOfElements), mean(0.), variance(0.),
				hostOutput(std::vector<spinor> (numberOfElements * testParameters.iterations)),	testParameters(testParameters){}

	~GaussianTester()
	{
		calculateMean();
		calculateVariance();

		kernelResult[0] = mean;
		kernelResult[1] = sqrt(variance);
	}

	double normalize(double valueIn)
	{
		return valueIn/= testParameters.iterations * numberOfElements * 24;
	}

	void calculateMean()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
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
protected:
	const int numberOfElements;
	double mean, variance;
	std::vector<spinor> hostOutput;
	const GaussianTestParameters & testParameters;
};

struct NonEvenOddGaussianSpinorfieldTestParameters : public GaussianTestParameters
{
	NonEvenOddGaussianSpinorfieldTestParameters(const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisonIn) :
		GaussianTestParameters(calculateReferenceValues_gaussian(), latticeExtendsIn, typeOfComparisonIn, false) {};
};

struct NonEvenGaussianSpinorfieldTester: public GaussianTester
{
	NonEvenGaussianSpinorfieldTester(const ParameterCollection parameterCollection, const GaussianTestParameters testParameters):
		GaussianTester("generate_gaussian_spinorfield", parameterCollection, testParameters, testParameters.getSpinorfieldSize() ){}
	~NonEvenGaussianSpinorfieldTester()
	{
		const hardware::buffers::Plain<spinor> outSpinor(spinorfieldElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			code->generate_gaussian_spinorfield_device(&outSpinor, prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct EvenOddGaussianSpinorfieldTestParameters : public GaussianTestParameters
{
	EvenOddGaussianSpinorfieldTestParameters(const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisonIn) :
		GaussianTestParameters(calculateReferenceValues_gaussian(), latticeExtendsIn, typeOfComparisonIn, true) {};
};

struct EvenOddGaussianSpinorfieldEvenOddTester: public GaussianTester
{
	EvenOddGaussianSpinorfieldEvenOddTester(const ParameterCollection parameterCollection, const GaussianTestParameters testParameters):
				GaussianTester("generate_gaussian_spinorfield_eo", parameterCollection, testParameters, testParameters.getEvenOddSpinorfieldSize() ) {};
	~EvenOddGaussianSpinorfieldEvenOddTester()
	{
		const hardware::buffers::Spinor outSpinor(spinorfieldEvenOddElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			code->generate_gaussian_spinorfield_eo_device(&outSpinor, prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct ConvertEvenOddTestParameters: public SpinorTestParameters
{
	const bool fillEvenSites;

	ConvertEvenOddTestParameters(LatticeExtents latticeExtendsIn, const bool fillEvenSitesIn):
		SpinorTestParameters(calculateReferenceValues_convert_eo(getEvenOddSpinorfieldSize(latticeExtendsIn), fillEvenSitesIn), latticeExtendsIn, true),
		fillEvenSites(fillEvenSitesIn) {}
};

struct ConvertToEvenOddTester: public SpinorTester
{
	ConvertToEvenOddTester(const ParameterCollection & parameterCollection, const ConvertEvenOddTestParameters testParameters):
	SpinorTester("convert_to_eo", parameterCollection, testParameters)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);

			in.load( createSpinorfieldWithOnesAndZerosDependingOnSiteParity( testParameters.fillEvenSites ) );
			code->convert_to_eoprec_device(&in2, &in3, &in) ;

			code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);
			code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
			doubleBuffer->dump(&kernelResult[1]);
		}
};

struct ConvertFromEvenOddTestParameters: public SpinorTestParameters
{
	const bool fillEvenSites;

	ConvertFromEvenOddTestParameters(LatticeExtents latticeExtendsIn, const bool fillEvenSitesIn):
		SpinorTestParameters( calculateReferenceValues_convertFromEvenOdd(getSpinorfieldSize(latticeExtendsIn) ), latticeExtendsIn, true),
		fillEvenSites(fillEvenSitesIn) {}
};

struct ConvertFromEvenOddTester: public SpinorTester
{
	ConvertFromEvenOddTester(const ParameterCollection & parameterCollection, const ConvertEvenOddTestParameters testParameters):
		SpinorTester("convert_from_eo", parameterCollection, testParameters)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);

			fillTwoSpinorBuffersDependingOnParity(&in2, &in3);
			code->convert_from_eoprec_device(&in2, &in3, &in);

			code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);
			}
};

struct SaxsbypzEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	SaxsbypzEvenOddTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_saxsbypz(getEvenOddSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 4){}
};

struct SaxsbypzEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxsbypzEvenOddTester(const ParameterCollection & parameterCollection, const SaxsbypzEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsResult("saxsbypz_eo", parameterCollection, testParameters)
		{
			code->saxsbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct SaxsbypzTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	SaxsbypzTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_saxsbypz(getSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 4){}
};

struct SaxsbypzTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxsbypzTester(const ParameterCollection & parameterCollection, const SaxsbypzTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("saxsbypz", parameterCollection, testParameters)
		{
			code->saxsbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct SquarenormTestParameters: public NonEvenOddLinearCombinationTestParameters
{
	SquarenormTestParameters(const LatticeExtents & latticeExtendsIn, const SpinorFillTypes & fillTypesIn) :
		NonEvenOddLinearCombinationTestParameters{calculateReferenceValues_globalSquarenorm( getSpinorfieldSize(latticeExtendsIn), fillTypesIn),
		latticeExtendsIn, fillTypesIn} {};
};

struct SquarenormTester: public NonEvenOddLinearCombinationTester3
{
	SquarenormTester(const ParameterCollection & parameterCollection, const SquarenormTestParameters & testParameters):
		NonEvenOddLinearCombinationTester3("global squarenorm", parameterCollection, testParameters) {}
};

struct SquarenormEvenOddTestParameters: public EvenOddLinearCombinationTestParameters
{
	SquarenormEvenOddTestParameters(const LatticeExtents & latticeExtendsIn, const SpinorFillTypes & fillTypesIn) :
		EvenOddLinearCombinationTestParameters{calculateReferenceValues_globalSquarenorm( getEvenOddSpinorfieldSize(latticeExtendsIn), fillTypesIn) , latticeExtendsIn, fillTypesIn} {};
};

struct SquarenormEvenOddTester: public EvenOddLinearCombinationTester3
{
	SquarenormEvenOddTester(const ParameterCollection & parameterCollection, const SquarenormEvenOddTestParameters & testParameters):
		EvenOddLinearCombinationTester3("global_squarenorm_eo", parameterCollection, testParameters) {}
};

struct ScalarProductTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	ScalarProductTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_scalarProduct( getSpinorfieldSize(latticeExtendsIn), fillTypesIn), latticeExtendsIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2){};
};

struct ScalarProductTester: public NonEvenOddLinearCombinationTester
{
	ScalarProductTester(const ParameterCollection & parameterCollection, const ScalarProductTestParameters testParameters):
		NonEvenOddLinearCombinationTester("scalar_product", parameterCollection, testParameters)
	{
		code->set_complex_to_scalar_product_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
		hmc_complex resultTmp;
		complexNums.at(0)->dump(&resultTmp);

		kernelResult[0] = resultTmp.re;
		kernelResult[1] = resultTmp.im;
	}
};

struct ScalarProductEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	ScalarProductEvenOddTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_scalarProduct(getEvenOddSpinorfieldSize(latticeExtendsIn), fillTypesIn), latticeExtendsIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2) {};
};

struct ScalarProductEvenOddTester: public EvenOddLinearCombinationTester
{
	ScalarProductEvenOddTester(const ParameterCollection & parameterCollection, const ScalarProductEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTester("scalar_product_eo", parameterCollection, testParameters)
	{
		code->set_complex_to_scalar_product_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
		hmc_complex resultTmp;
		complexNums.at(0)->dump(&resultTmp);

		kernelResult.at(0) = resultTmp.re;
		kernelResult.at(1) = resultTmp.im;
	}
};

struct ZeroTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	ZeroTestParameters(LatticeExtents latticeExtendsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_zero(), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

struct ZeroTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	ZeroTester(const ParameterCollection & parameterCollection, const ZeroTestParameters & testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("zero", parameterCollection, testParameters)
		{
			code->set_zero_spinorfield_device(getOutSpinor());
		}
};


struct ColdTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	ColdTestParameters(LatticeExtents latticeExtendsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_cold(false), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

struct ColdTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	ColdTester(const ParameterCollection & parameterCollection, const ColdTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("cold", parameterCollection, testParameters)
		{
			code->set_spinorfield_cold_device(getOutSpinor());
		}
};

struct ZeroEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	ZeroEvenOddTestParameters(LatticeExtents latticeExtendsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_zero(), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

struct ZeroEvenOddTester: public EvenOddLinearCombinationTester
{
	ZeroEvenOddTester(const ParameterCollection & parameterCollection, const ZeroEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTester ("zero_eo", parameterCollection, testParameters)
		{
			code->set_zero_spinorfield_eoprec_device(getOutSpinor());
		}
};

struct ColdEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	ColdEvenOddTestParameters(LatticeExtents latticeExtendsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_cold(true), latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}) {};
};

struct ColdEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
{
	ColdEvenOddTester(const ParameterCollection & parameterCollection, const ColdEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsResult("cold_eo",parameterCollection, testParameters)
		{
			code->set_eoprec_spinorfield_cold_device(getOutSpinor());
		}
};

struct SaxTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	SaxTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_sax(getSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 2){}
};

struct SaxTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxTester(const ParameterCollection & parameterCollection, const SaxTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("sax", parameterCollection, testParameters)
		{
			code->sax_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	SaxEvenOddTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_sax(getEvenOddSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 2){}
};

struct SaxEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxEvenOddTester(const ParameterCollection & parameterCollection, const SaxEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsResult("sax_eo", parameterCollection, testParameters)
		{
			code->sax_eoprec_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpyTestParameters : public NonEvenOddLinearCombinationTestParameters
{
	SaxpyTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		NonEvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpy(getSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 3){}
};

struct SaxpyTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxpyTester(const ParameterCollection & parameterCollection, const SaxpyTestParameters testParameters):
	NonEvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy", parameterCollection, testParameters)
		{
			code->saxpy_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpyArgTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxpyArgTester(const ParameterCollection & parameterCollection, const SaxpyTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_arg", parameterCollection, testParameters)
		{
			code->saxpy_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), getOutSpinor());
		}
};
struct SaxpyEvenOddTestParameters : public EvenOddLinearCombinationTestParameters
{
	SaxpyEvenOddTestParameters(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
		EvenOddLinearCombinationTestParameters(calculateReferenceValues_saxpy(getEvenOddSpinorfieldSize(latticeExtendsIn), coefficientsIn), latticeExtendsIn, fillTypesIn, coefficientsIn, 3){}
};

struct SaxpyEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxpyEvenOddTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_eo", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpyArgEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsResult
{
	SaxpyArgEvenOddTester(const ParameterCollection & parameterCollection, const SaxpyEvenOddTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsResult("saxpy_arg_eo", parameterCollection, testParameters)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(0), testParameters.coefficients.at(0), getOutSpinor());
		}
};

void testEvenOddScalarProduct( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (lE, sF);
}

void testNonEvenOddSax(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxTester, SaxTestParameters> (lE, cN);
}

void testEvenOddSax(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (lE, cN);
}

void testNonEvenOddSaxpy(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyTester, SaxpyTestParameters> (lE, cN);
}

void testNonEvenOddSaxpyArg(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyArgTester, SaxpyTestParameters> (lE, cN);
}

void testEvenOddSaxpy(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (lE, cN);
}

void testEvenOddSaxpyArg(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (lE, cN);
}

void testEvenOddGaussianSpinorfield( const LatticeExtents lE)
{
	performTest<NonEvenGaussianSpinorfieldTester, NonEvenOddGaussianSpinorfieldTestParameters>(lE, ComparisonType::smallerThan);
}

void testNonEvenOddGaussianSpinorfield( const LatticeExtents lE)
{
	performTest<EvenOddGaussianSpinorfieldEvenOddTester, EvenOddGaussianSpinorfieldTestParameters>(lE, ComparisonType::smallerThan);
}

void testConvertToEvenOdd(const LatticeExtents lE, const bool fillEvenSites)
{
	performTest<ConvertToEvenOddTester, ConvertEvenOddTestParameters> (lE, fillEvenSites);
}

void testConvertFromEvenOdd(const LatticeExtents lE, const bool fillEvenSites)
{
	performTest<ConvertFromEvenOddTester, ConvertEvenOddTestParameters> (lE, fillEvenSites);
}

void testSaxsbypzEvenOdd( const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (lE, cN);
}
void testNonEvenOddSaxsbypz(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxsbypzTester, SaxsbypzTestParameters> (lE, cN);
}

void testNonEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<SquarenormTester, SquarenormTestParameters> (lE, sF);
}

void testEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (lE, sF);
}

void testNonEvenOddScalarProduct( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ScalarProductTester, ScalarProductTestParameters> (lE, sF);
}


BOOST_AUTO_TEST_SUITE(SPINORTESTER_BUILD)

	BOOST_AUTO_TEST_CASE( BUILDFROMPARAMETERS )
	{
		const hardware::HardwareParametersMockup hardwareParameters(4,4);
		const hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(4,4);
		ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
		const SpinorTestParameters testParameters;
		BOOST_CHECK_NO_THROW( SpinorTester( "build all kernels", parameterCollection, testParameters) );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GLOBAL_SQUARENORM)

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_1 )
	{
		testNonEvenOddSquarenorm(  LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_2 )
	{
		testNonEvenOddSquarenorm( LatticeExtents{ns8, nt4}, SpinorFillTypes{ SpinorFillType::ascendingComplex} );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_1 )
	{
		testNonEvenOddSquarenorm( LatticeExtents{ns8, nt12}, SpinorFillTypes{ SpinorFillType::one } );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_2 )
	{
		testNonEvenOddSquarenorm( LatticeExtents{ns12, nt16}, SpinorFillTypes{ SpinorFillType::one } );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_3 )
	{
		testNonEvenOddSquarenorm( LatticeExtents{ns16, nt8}, SpinorFillTypes{ SpinorFillType::one} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GLOBAL_SQUARENORM_EO)

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns16, nt8}, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns16, nt8}, SpinorFillTypes{ SpinorFillType::ascendingComplex });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns4, nt16}, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns8, nt4}, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns16, nt16}, SpinorFillTypes { SpinorFillType::one });
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_3 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_4 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns8, nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns8, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns8, nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns8, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns4, nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_3 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns8, nt8}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_4 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns16, nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns16, nt8}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		testEvenOddScalarProduct(LatticeExtents{ns4, nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ZERO)

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		performTest<ZeroTester, ZeroTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD)

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		performTest<ColdTester, ColdTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ZERO_EO)

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		performTest<ZeroEvenOddTester, ZeroEvenOddTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_EO)

	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		performTest<ColdEvenOddTester, ColdEvenOddTestParameters> (LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX)

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		testNonEvenOddSax (LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		testNonEvenOddSax(LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		testNonEvenOddSax(LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		testNonEvenOddSax(LatticeExtents{ns16, nt8}, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	BOOST_AUTO_TEST_CASE( SAX_EO_1 )
	{
		testEvenOddSax(LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_2 )
	{
		testEvenOddSax(LatticeExtents{ns4, nt8}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_3 )
	{
		testEvenOddSax(LatticeExtents{ns8, nt4}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_4 )
	{
		testEvenOddSax(LatticeExtents{ns16, nt8}, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		testNonEvenOddSaxpy(LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		testNonEvenOddSaxpy(LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		testNonEvenOddSaxpy(LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		testNonEvenOddSaxpy(LatticeExtents{ns8, nt8}, ComplexNumbers {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
		testNonEvenOddSaxpyArg(LatticeExtents{ns12, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
		testNonEvenOddSaxpyArg(LatticeExtents{ns4, nt12}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
		testNonEvenOddSaxpyArg(LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
		testNonEvenOddSaxpyArg(LatticeExtents{ns8, nt16}, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)


	BOOST_AUTO_TEST_CASE( SAXPY_EO_1 )
	{
		testEvenOddSaxpy(LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_2 )
	{
		testEvenOddSaxpy(LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_3 )
	{
		testEvenOddSaxpy(LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_4 )
	{
		testEvenOddSaxpy(LatticeExtents{ns8, nt8}, ComplexNumbers {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_5 )
	{
		testEvenOddSaxpyArg(LatticeExtents{ns12, nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_6 )
	{
		testEvenOddSaxpyArg(LatticeExtents{ns4, nt12}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_7 )
	{
		testEvenOddSaxpyArg(LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_8 )
	{
		testEvenOddSaxpyArg(LatticeExtents{ns8, nt16}, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ)

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_1 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		testNonEvenOddSaxsbypz(LatticeExtents{ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
	{
		testSaxsbypzEvenOdd(LatticeExtents{ns8, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
		testConvertToEvenOdd(LatticeExtents{ns4, nt4}, true);
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
		testConvertToEvenOdd(LatticeExtents{ns8, nt4}, false);
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		testConvertFromEvenOdd(LatticeExtents{ns4, nt4}, true);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_NONEO, 1)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_NONEO )
	{
		testNonEvenOddGaussianSpinorfield( LatticeExtents{ns8, nt4} );
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_EO, 1)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO )
	{
		testEvenOddGaussianSpinorfield(LatticeExtents{ns12, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()
