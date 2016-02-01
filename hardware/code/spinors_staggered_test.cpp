/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra,
 * Francesca Cuteri
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
#include "PrngSpinorTester.hpp"

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

struct LinearCombinationTestParameters : public SpinorStaggeredTestParameters
{
	LinearCombinationTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn) :
		TestParameters(latticeExtentsIn), SpinorStaggeredTestParameters(latticeExtentsIn, fillTypesIn), coefficients(ComplexNumbers{{1.,0.}}), numberOfSpinors(1) {};
	LinearCombinationTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers cN, const size_t numberOfSpinorsIn) :
		TestParameters(latticeExtentsIn), SpinorStaggeredTestParameters(latticeExtentsIn, fillTypesIn), coefficients(cN), numberOfSpinors(numberOfSpinorsIn) {};
	const ComplexNumbers coefficients;
	const NumberOfSpinors numberOfSpinors;
};

struct SaxVecAndSqnormEvenOddTestParameters: public SpinorStaggeredTestParameters
{
	SaxVecAndSqnormEvenOddTestParameters(const LatticeExtents latticeExtentsIn, const ComplexNumbers coefficientsIn, const int numEqsIn):
		TestParameters(latticeExtentsIn), SpinorStaggeredTestParameters(latticeExtentsIn),
		coefficients(coefficientsIn),
		numEqs(numEqsIn){};
		ComplexNumbers coefficients;
		const int numEqs;
};

struct GaussianTestParameters: public SpinorStaggeredTestParameters
{
	GaussianTestParameters(const LatticeExtents latticeExtentsIn, const ComparisonType & typeOfComparisonIn) :
		TestParameters(latticeExtentsIn, ComparisonType::smallerThan), SpinorStaggeredTestParameters(latticeExtentsIn, typeOfComparisonIn),	iterations(1000){};

	const unsigned int iterations;
};

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest, const bool needEvenOdd)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, needEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const bool needEvenOdd)
{
	LinearCombinationTestParameters parametersForThisTest(latticeExtendsIn, fillTypesIn);
	callTest<TesterClass, LinearCombinationTestParameters>(parametersForThisTest, needEvenOdd);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const ComplexNumbers alphaIn, const int numberOfSpinors, const bool needEvenOdd )
{
	LinearCombinationTestParameters parametersForThisTest(latticeExtendsIn, SpinorFillTypes{SpinorFillType::ascendingComplex}, alphaIn, numberOfSpinors);
	callTest<TesterClass, LinearCombinationTestParameters>(parametersForThisTest, needEvenOdd);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const ComplexNumbers alphaIn, const int numEqsIn )
{
	SaxVecAndSqnormEvenOddTestParameters parametersForThisTest(latticeExtendsIn, alphaIn, numEqsIn);
	callTest<TesterClass, SaxVecAndSqnormEvenOddTestParameters>(parametersForThisTest, true);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const SpinorFillTypes sF, const ComplexNumbers alphaIn, const int numberOfSpinors, const bool needEvenOdd )
{
	LinearCombinationTestParameters parametersForThisTest(latticeExtendsIn, sF, alphaIn, numberOfSpinors);
	callTest<TesterClass, LinearCombinationTestParameters>(parametersForThisTest, needEvenOdd);
}

template<typename TesterClass> void performTest(const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisonIn, const bool needEvenOdd )
{
	GaussianTestParameters parametersForThisTest(latticeExtendsIn, typeOfComparisonIn);
	callTest<TesterClass, GaussianTestParameters>(parametersForThisTest, needEvenOdd);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const bool fillEvenSitesIn )
{
	SpinorStaggeredTestParameters parametersForThisTest(latticeExtendsIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest, fillEvenSitesIn);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn)
{
	SpinorStaggeredTestParameters parametersForThisTest(latticeExtendsIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

struct NonEvenOddLinearCombinationTester : public NonEvenOddSpinorStaggeredTester
{
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		NonEvenOddSpinorStaggeredTester(kernelName, pC, tP, rV(calculateSpinorfieldSize(tP.latticeExtents), tP.coefficients))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		NonEvenOddSpinorStaggeredTester(kernelName, pC, tP, rV)
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
protected:
	std::vector<const hardware::buffers::Plain<su3vec> *> spinorfields;
	std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
	const hardware::buffers::Plain<su3vec> * getOutSpinor() const
	{
		return spinorfields.back();
	}
private:
	void loadCoefficients(const LinearCombinationTestParameters & testParameters)
	{
		for (auto coefficient : testParameters.coefficients)
		{
			complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
			complexNums.back()->load(&coefficient);
		}
	}
	void loadSpinorfields(const LinearCombinationTestParameters & tP)
	{
		for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
		{
			NonEvenOddSpinorStaggeredfieldCreator ssf(tP.latticeExtents);
			spinorfields.push_back(new hardware::buffers::Plain<su3vec>(calculateSpinorfieldSize(tP.latticeExtents), device));
			(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(ssf.createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(ssf.createSpinorfield(tP.fillTypes.at(number)));
		}
	}
};

struct NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult : public NonEvenOddLinearCombinationTester
{
	NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const SpinorFillTypes)):
		NonEvenOddLinearCombinationTester(kernelName, pC, tP, rV(calculateSpinorfieldSize(tP.latticeExtents), tP.fillTypes) ) {}
	NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		NonEvenOddLinearCombinationTester(kernelName, pC, tP, rV) {}
	NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		NonEvenOddLinearCombinationTester(kernelName, pC, tP, rV) {}
	~NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult()
	{
		calcSquarenormAndStoreAsKernelResult(getOutSpinor());
	}
};

struct EvenOddLinearCombinationTester : public EvenOddSpinorStaggeredTester
{
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const SpinorFillTypes)):
		EvenOddSpinorStaggeredTester(kernelName, pC, tP, rV( calculateEvenOddSpinorfieldSize(tP.latticeExtents), tP.fillTypes))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		EvenOddSpinorStaggeredTester(kernelName, pC, tP, rV( calculateEvenOddSpinorfieldSize(tP.latticeExtents), tP.coefficients))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		EvenOddSpinorStaggeredTester(kernelName, pC, tP, rV)
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	protected:
		std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
		std::vector<const hardware::buffers::SU3vec *> spinorfields;
		const hardware::buffers::SU3vec * getOutSpinor() const
		{
			return spinorfields.back();
		}
	private:
		void loadCoefficients(const LinearCombinationTestParameters & testParameters)
		{
			for (auto coefficient : testParameters.coefficients)
			{
				complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
				complexNums.back()->load(&coefficient);
			}
		}
		void loadSpinorfields(const LinearCombinationTestParameters & tP)
		{
			for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
			{
				EvenOddSpinorStaggeredfieldCreator ssf(tP.latticeExtents);
				spinorfields.push_back(new hardware::buffers::SU3vec(tP.latticeExtents, device));
				(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(ssf.createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(ssf.createSpinorfield(tP.fillTypes.at(number)));
			}
		}
};

struct EvenOddLinearCombinationTesterWithSquarenormAsKernelResult : public EvenOddLinearCombinationTester
{
	EvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const SpinorFillTypes)):
		EvenOddLinearCombinationTester(kernelName, pC, tP, rV) {}
	EvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		EvenOddLinearCombinationTester(kernelName, pC, tP, rV) {}
	EvenOddLinearCombinationTesterWithSquarenormAsKernelResult(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		EvenOddLinearCombinationTester(kernelName, pC, tP, rV) {}
	~EvenOddLinearCombinationTesterWithSquarenormAsKernelResult()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinor());
	}
};

struct SquarenormTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SquarenormTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters & testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("squarenorm", parameterCollection, testParameters, calculateReferenceValues_squarenorm) {}
};

struct ScalarProductTester: public NonEvenOddLinearCombinationTester
{
	ScalarProductTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTester("scalar_product", parameterCollection, testParameters, calculateReferenceValues_scalarProduct(calculateSpinorfieldSize(testParameters.latticeExtents), testParameters.fillTypes))
	{
		code->set_complex_to_scalar_product_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
		hmc_complex resultHost;
		complexNums.at(0)->dump(&resultHost);

		kernelResult[0] = resultHost.re;
		kernelResult[1] = resultHost.im;
	}
};

struct ZeroTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ZeroTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters & testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("zero", parameterCollection, testParameters, calculateReferenceValues_zero())
		{
			code->set_zero_spinorfield_device(getOutSpinor());
		}
};

struct ColdTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ColdTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("cold", parameterCollection, testParameters, calculateReferenceValues_cold(false))
		{
			code->set_cold_spinorfield_device(getOutSpinor());
		}
};

struct SaxTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			code->sax_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpyTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpbypzTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbypzTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpbypz", parameterCollection, testParameters, calculateReferenceValues_saxpbypz)
		{
			code->saxpbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(3), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct ConvertToEvenOddTester: public SpinorStaggeredTester2
{
	ConvertToEvenOddTester(const ParameterCollection & parameterCollection, SpinorStaggeredTestParameters testParameters, const bool fillEvenSitesIn):
		SpinorStaggeredTester2("convert_to_eo", parameterCollection, testParameters,
				calculateReferenceValues_convert_eo(calculateEvenOddSpinorfieldSize(testParameters.latticeExtents), fillEvenSitesIn))
		{
			const hardware::buffers::Plain<su3vec> in(calculateSpinorfieldSize(testParameters.latticeExtents), device);
			const hardware::buffers::SU3vec in2(testParameters.latticeExtents, device);
			const hardware::buffers::SU3vec in3(testParameters.latticeExtents, device);

			NonEvenOddSpinorStaggeredfieldCreator ssf(testParameters.latticeExtents);
			in.load(ssf.createSpinorfieldWithOnesAndZerosDependingOnSiteParity( fillEvenSitesIn ));
			code->convert_to_eoprec_device(&in2, &in3, &in);

			code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);
			code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
			doubleBuffer->dump(&kernelResult[1]);
		}
};

struct ConvertFromEvenOddTester: public SpinorStaggeredTester2
{
	ConvertFromEvenOddTester(const ParameterCollection & parameterCollection, SpinorStaggeredTestParameters testParameters, const bool fillEvenSitesIn):
		SpinorStaggeredTester2("convert_from_eo", parameterCollection, testParameters,
				calculateReferenceValues_convertFromEvenOdd(calculateSpinorfieldSize(testParameters.latticeExtents)))
		{
			const hardware::buffers::Plain<su3vec> in(calculateSpinorfieldSize(testParameters.latticeExtents), device);
			const hardware::buffers::SU3vec in2(testParameters.latticeExtents, device);
			const hardware::buffers::SU3vec in3(testParameters.latticeExtents, device);

			EvenOddSpinorStaggeredfieldCreator ssf(testParameters.latticeExtents);
			in2.load(ssf.createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity( fillEvenSitesIn ));
			in3.load(ssf.createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity( fillEvenSitesIn ));
			code->convert_from_eoprec_device(&in2, &in3, &in);
			code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);

		}
};

struct SquarenormEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SquarenormEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters & testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("squarenorm_eo", parameterCollection, testParameters, calculateReferenceValues_squarenorm){}
};

struct ScalarProductEvenOddRealTester: public EvenOddLinearCombinationTester
{
	ScalarProductEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTester("scalar_product_real_eo", parameterCollection, testParameters, calculateReferenceValues_scalarProduct)
		{
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);
			code->set_float_to_scalar_product_real_part_eoprec_device(spinorfields.at(0), spinorfields.at(1), &sqnorm);
			hmc_float resultTmp;
			sqnorm.dump(&resultTmp);
			kernelResult[0] = resultTmp;
			kernelResult[1] = 0;
		}
};

struct ScalarProductEvenOddComplexTester: public EvenOddLinearCombinationTester
{
	ScalarProductEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTester("scalar_product_real_eo", parameterCollection, testParameters, calculateReferenceValues_scalarProduct)
		{
			code->set_complex_to_scalar_product_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
			hmc_complex resultTmp;
			complexNums.at(0)->dump(&resultTmp);
			kernelResult[0] = resultTmp.re;
			kernelResult[1] = resultTmp.im;
		}
};

struct ZeroEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ZeroEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("zero_eo", parameterCollection, testParameters, calculateReferenceValues_zero())
		{
			code->set_zero_spinorfield_eoprec_device(getOutSpinor());
		}
};

struct ColdEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ColdEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("cold_eo", parameterCollection, testParameters, calculateReferenceValues_cold(true))
		{
			code->set_cold_spinorfield_eoprec_device(getOutSpinor());
		}
};

struct SaxEvenOddComplexTester: EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_eo_complex", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			code->sax_eoprec_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxArgEvenOddComplexTester: EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_arg_eo_complex", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			code->sax_eoprec_device(spinorfields.at(0), testParameters.coefficients.at(0), getOutSpinor());
		}
};

struct SaxEvenOddRealTester: EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_eo_real", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			hmc_complex complexNumbers;
			complexNums.at(0)->dump(&complexNumbers);
			code->sax_eoprec_device(spinorfields.at(0), complexNumbers.re, getOutSpinor());
		}
};

struct SaxArgEvenOddRealTester: EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxArgEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_arg_eo_real", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			code->sax_eoprec_device(spinorfields.at(0), testParameters.coefficients.at(0).re, getOutSpinor());
		}
};

struct SaxVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxVecEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_eo_vec", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
			std::vector<hmc_float> alpha_host_real_vec(5, testParameters.coefficients.at(0).re);
			const int index_alpha = 3;
			alpha_real_vec.load(&alpha_host_real_vec[0]);
			code->sax_eoprec_device(spinorfields.at(0), &alpha_real_vec, index_alpha, getOutSpinor());
		}
};

struct SaxpyEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_eo_complex", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), getOutSpinor());
		}
};

struct SaxpyArgEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_arg_eo_complex", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), getOutSpinor());
		}
};

struct SaxpyEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_eo_real", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			hmc_complex complexNumbers;
			complexNums.at(0)->dump(&complexNumbers);
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNumbers.re, getOutSpinor());
		}
};

struct SaxpyArgEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyArgEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_arg_eo_real", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0).re, getOutSpinor());
		}
};

struct SaxpyVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyVecEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_eo_vec", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
			std::vector<hmc_float> alpha_host_real_vec(5, testParameters.coefficients.at(0).re);
			const int index_alpha = 3;
			alpha_real_vec.load(&alpha_host_real_vec[0]);
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(1), &alpha_real_vec, index_alpha, getOutSpinor());
		}
};

struct SaxpbyEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbyEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpby_eo_complex", parameterCollection, testParameters, calculateReferenceValue_saxpby)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), testParameters.coefficients.at(1), getOutSpinor());
		}
};

struct SaxpbyArgEvenOddComplexTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbyArgEvenOddComplexTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpby_arg_eo_complex", parameterCollection, testParameters, calculateReferenceValue_saxpby)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0), testParameters.coefficients.at(1), getOutSpinor());
		}
};

struct SaxpbyEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbyEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpby_eo_real", parameterCollection, testParameters, calculateReferenceValue_saxpby)
		{
			hmc_complex complexNumbers1;
			hmc_complex complexNumbers2;
			complexNums.at(0)->dump(&complexNumbers1);
			complexNums.at(1)->dump(&complexNumbers2);
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNumbers1.re, complexNumbers2.re, getOutSpinor());
		}
};

struct SaxpbyArgEvenOddRealTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbyArgEvenOddRealTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpby_arg_eo_real", parameterCollection, testParameters, calculateReferenceValue_saxpby)
		{
		    code->saxpby_eoprec_device(spinorfields.at(0), spinorfields.at(1), testParameters.coefficients.at(0).re, testParameters.coefficients.at(1).re, getOutSpinor());
		}
};

struct SaxpbyVecEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbyVecEvenOddTester(const ParameterCollection parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpby_eo_vec", parameterCollection, testParameters, calculateReferenceValue_saxpby)
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

struct SaxpbypzEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbypzEvenOddTester(const ParameterCollection & parameterCollection, LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpbypz_eo", parameterCollection, testParameters, calculateReferenceValues_saxpbypz)
		{
			code->saxpbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct SaxpbypzArgEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpbypzArgEvenOddTester(const ParameterCollection & parameterCollection, LinearCombinationTestParameters testParameters):
		 EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpbypz_arg_eo", parameterCollection, testParameters, calculateReferenceValues_saxpbypz)
		{
			code->saxpbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), testParameters.coefficients.at(0), testParameters.coefficients.at(1), getOutSpinor());
		}
};

struct SaxVecAndSqnormEvenOddTester: public EvenOddSpinorStaggeredTester
{
	SaxVecAndSqnormEvenOddTester(const ParameterCollection & parameterCollection, const SaxVecAndSqnormEvenOddTestParameters testParameters):
		EvenOddSpinorStaggeredTester("sax_vectorized_and_squarenorm_eoprec", parameterCollection, testParameters, calculateReferenceValue_sax_vec_and_sqnorm(calculateEvenOddSpinorfieldSize(testParameters.latticeExtents), testParameters.coefficients, testParameters.numEqs))
		{
			int NUM_EQS = testParameters.numEqs;
			EvenOddSpinorStaggeredfieldCreator ssf(testParameters.latticeExtents);
			const hardware::buffers::SU3vec in(testParameters.latticeExtents, device, NUM_EQS);
			const hardware::buffers::SU3vec out(testParameters.latticeExtents, device, NUM_EQS);
			in.load(ssf.createSpinorfield(SpinorFillType::one));

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

struct PrngTester: public SpinorStaggeredTester2
{
	PrngTester(const std::string kernelName, const ParameterCollection parameterCollection, const GaussianTestParameters & testParameters, const ReferenceValues & referenceValues):
		SpinorStaggeredTester2(kernelName, parameterCollection, testParameters, referenceValues),
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
		PrngTester(kernelName, parameterCollection, testParameters, calculateReferenceValues_gaussian()),
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
				if(i%100==0) logger.info() << "Run kernel for the " << i << "th time";
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

struct NonEvenOddGaussianStaggeredSpinorfieldTester: public GaussianTester
{
	NonEvenOddGaussianStaggeredSpinorfieldTester(const ParameterCollection & parameterCollection, const GaussianTestParameters testParameters):
		GaussianTester("generate_gaussian_staggeredspinorfield", parameterCollection, testParameters, calculateSpinorfieldSize(testParameters.latticeExtents)) {}
	~NonEvenOddGaussianStaggeredSpinorfieldTester()
	{
		const hardware::buffers::Plain<su3vec> outSpinor(numberOfElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++){
			code->set_gaussian_spinorfield_device(&outSpinor,prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct EvenOddGaussianStaggeredSpinorfieldTester: public GaussianTester
{
	EvenOddGaussianStaggeredSpinorfieldTester(const ParameterCollection & parameterCollection, const GaussianTestParameters testParameters):
		GaussianTester("generate_gaussian_staggeredspinorfield", parameterCollection, testParameters, calculateEvenOddSpinorfieldSize(testParameters.latticeExtents)) {}
	~EvenOddGaussianStaggeredSpinorfieldTester()
	{
		const hardware::buffers::SU3vec outSpinor(numberOfElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++){
			code->set_gaussian_spinorfield_eoprec_device(&outSpinor, prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

void testNonEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF )
{
	performTest<SquarenormTester> ( lE, sF, false );
}

void testNonEvenOddScalarProduct( const LatticeExtents lE, const SpinorFillTypes sF )
{
	performTest<ScalarProductTester> ( lE, sF, ComplexNumbers{{1,0}}, 2, false );
}

void testNonEvenOddZero( const LatticeExtents lE )
{
	performTest<ZeroTester> ( lE, ComplexNumbers{{1,0}}, 2, false );
}

void testNonEvenOddCold( const LatticeExtents lE )
{
	performTest<ColdTester> ( lE, ComplexNumbers{{1,0}}, 2, false );
}

void testNonEvenOddSax( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxTester> (lE, cN, 2, false);
}

void testNonEvenOddSaxpy( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyTester> ( lE, cN, 2, false );
}

void testNonEvenOddSaxpbypz( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbypzTester> ( lE, cN, 4, false );
}

void testConvertToEvenOdd( const LatticeExtents lE, const bool fillEvenSites)
{
	performTest<ConvertToEvenOddTester> (lE, fillEvenSites);
}

void testConvertFromEvenOdd( const LatticeExtents lE, const bool fillEvenSites)
{
	performTest<ConvertFromEvenOddTester> (lE, fillEvenSites);
}

void testNonEvenOddGaussianSpinorfield( const LatticeExtents lE)
{
	performTest<NonEvenOddGaussianStaggeredSpinorfieldTester>(lE, ComparisonType::smallerThan, false);
}

void testEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF )
{
	performTest<SquarenormEvenOddTester> ( lE, sF, true );
}

void testEvenOddScalarProductComplex( const LatticeExtents lE, const SpinorFillTypes sF )
{
	performTest<ScalarProductEvenOddComplexTester> ( lE, sF, ComplexNumbers{{1,0}}, 2, true );
}

void testEvenOddScalarProductReal( const LatticeExtents lE, const SpinorFillTypes sF )
{
	performTest<ScalarProductEvenOddRealTester> ( lE, sF, ComplexNumbers{{1,0}}, 2, true );
}

void testEvenOddZero( const LatticeExtents lE )
{
	performTest<ZeroEvenOddTester> ( lE, ComplexNumbers{{1,0}}, 2, true );
}

void testEvenOddCold( const LatticeExtents lE )
{
	performTest<ColdEvenOddTester> ( lE, ComplexNumbers{{1,0}}, 2, true );
}

void testEvenOddSaxComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxEvenOddComplexTester> ( lE, cN, 2, true );
}

void testEvenOddSaxArgComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxArgEvenOddComplexTester> ( lE, cN, 2, true );
}

void testEvenOddSaxReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxEvenOddRealTester> ( lE, cN, 2, true );
}

void testEvenOddSaxArgReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxArgEvenOddRealTester> ( lE, cN, 2, true );
}

void testEvenOddSaxVecReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxVecEvenOddTester> ( lE, cN, 2, true );
}

void testEvenOddSaxpyComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyEvenOddComplexTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpyArgComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyArgEvenOddComplexTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpyReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyEvenOddRealTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpyArgReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyArgEvenOddRealTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpyVecReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpyVecEvenOddTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbyComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyEvenOddComplexTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbyArgComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyArgEvenOddComplexTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbyReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyEvenOddRealTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbyArgReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyArgEvenOddRealTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbyVecReal( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyVecEvenOddTester> ( lE, cN, 3, true );
}

void testEvenOddSaxpbypzComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbypzEvenOddTester> ( lE, cN, 4, true );
}

void testEvenOddSaxpbypzArgComplex( const LatticeExtents lE, const ComplexNumbers cN )
{
	performTest<SaxpbyArgEvenOddComplexTester> ( lE, cN, 4, true );
}

void testEvenOddGaussianSpinorfield( const LatticeExtents lE)
{
	performTest<EvenOddGaussianStaggeredSpinorfieldTester>(lE, ComparisonType::smallerThan, true);
}

void testEvenOddSaxVecAndSqnorm( const LatticeExtents lE, const ComplexNumbers cN, const int numEqs)
{
	performTest<SaxVecAndSqnormEvenOddTester> (lE, cN, numEqs);
}

BOOST_AUTO_TEST_SUITE(SQUARENORM)

	BOOST_AUTO_TEST_CASE( SQUARENORM_1 )
	{
	   testNonEvenOddSquarenorm(  LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_2 )
	{
	   testNonEvenOddSquarenorm(  LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::ascendingComplex} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_1 )
	{
		testNonEvenOddSquarenorm(  LatticeExtents{ns8, nt8}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_2 )
	{
		testNonEvenOddSquarenorm(  LatticeExtents{ns12, nt12}, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_3 )
	{
		testNonEvenOddSquarenorm(  LatticeExtents{ns16, nt16}, SpinorFillTypes{ SpinorFillType::one} );
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
	testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns4, nt4}, SpinorFillTypes{ SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex} );
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns8, nt8}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns12, nt12}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		testNonEvenOddScalarProduct( LatticeExtents{ns16, nt16}, SpinorFillTypes{ SpinorFillType::one, SpinorFillType::one} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ZERO)

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		testNonEvenOddZero( LatticeExtents{ns4, nt4} );
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(COLD)

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		testNonEvenOddCold( LatticeExtents{ns4, nt4} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX)

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		testNonEvenOddSax( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}} );
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		testNonEvenOddSax( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}} );
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		testNonEvenOddSax( LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}} );
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		testNonEvenOddSax( LatticeExtents{ns16,nt4}, ComplexNumbers {{1.,1.}} );
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		testNonEvenOddSaxpy( LatticeExtents{ns4, nt4}, ComplexNumbers{{0.,0.}} );
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		testNonEvenOddSaxpy( LatticeExtents{ns8,nt4}, ComplexNumbers{{1.,0.}} );
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		testNonEvenOddSaxpy (LatticeExtents{ns4,nt8}, ComplexNumbers{{0.,1.}} );
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		testNonEvenOddSaxpy ( LatticeExtents{ns8,nt8}, ComplexNumbers{{1.,1.}} );

	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPBYPZ)

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_1 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns4,nt4}, ComplexNumbers{{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns8, nt16}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		testNonEvenOddSaxpbypz( LatticeExtents{ns16, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
	testNonEvenOddGaussianSpinorfield( LatticeExtents{ns8,nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
	    testConvertToEvenOdd( LatticeExtents {ns4, nt4}, true );
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
		testConvertToEvenOdd( LatticeExtents {ns8, nt4}, false );
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		testConvertFromEvenOdd( LatticeExtents{ns4, nt4}, true );
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
		testConvertFromEvenOdd( LatticeExtents{ns8, nt4}, false );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SQUARENORM_EO)

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		testEvenOddSquarenorm( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns8,nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns4,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_3 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_4 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		testEvenOddScalarProductComplex  ( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_1 )
	{
		testEvenOddScalarProductReal  ( LatticeExtents{ns8,nt4}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_2 )
	{
		testEvenOddScalarProductReal  ( LatticeExtents{ns4,nt4}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_1 )
	{
		testEvenOddScalarProductReal  ( LatticeExtents{ns8,nt8}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_2 )
	{
		testEvenOddScalarProductReal  ( LatticeExtents{ns12,nt12}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_REDUCTION_3 )
	{
		testEvenOddScalarProductReal  ( LatticeExtents{ns16,nt16}, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ZERO_EO)

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		testEvenOddZero( LatticeExtents {ns4, nt4} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_EO)

	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		testEvenOddCold( LatticeExtents{ns4,nt4} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_1 )
	{
	    testEvenOddSaxComplex( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_2 )
	{
	    testEvenOddSaxComplex( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_3 )
	{
	    testEvenOddSaxComplex( LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_4 )
	{
	    testEvenOddSaxComplex( LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_1 )
	{
	    testEvenOddSaxArgComplex( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_2 )
	{
	    testEvenOddSaxArgComplex( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_3 )
	{
	    testEvenOddSaxArgComplex( LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_CPLX_EO_4 )
	{
	    testEvenOddSaxArgComplex( LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_1 )
	{
	    testEvenOddSaxReal( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_2 )
	{
	    testEvenOddSaxReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_3 )
	{
	    testEvenOddSaxReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_1 )
	{
	    testEvenOddSaxArgReal( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_2 )
	{
	    testEvenOddSaxArgReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_ARG_REAL_EO_3 )
	{
	    testEvenOddSaxArgReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_1 )
	{
	    testEvenOddSaxVecReal( LatticeExtents{ns4,nt4}, ComplexNumbers{{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_2 )
	{
	    testEvenOddSaxVecReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_3 )
	{
	    testEvenOddSaxVecReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_1 )
	{
	    testEvenOddSaxpyComplex( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_2 )
	{
	    testEvenOddSaxpyComplex( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_3 )
	{
	    testEvenOddSaxpyComplex( LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_4 )
	{
	    testEvenOddSaxpyComplex( LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_1 )
	{
	    testEvenOddSaxpyArgComplex( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_2 )
	{
	    testEvenOddSaxpyArgComplex( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_3 )
	{
	    testEvenOddSaxpyArgComplex( LatticeExtents{ns4,nt8}, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_CPLX_EO_4 )
	{
	    testEvenOddSaxpyArgComplex( LatticeExtents{ns16,nt8}, ComplexNumbers {{1.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_1 )
	{
	    testEvenOddSaxpyReal( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_2 )
	{
	    testEvenOddSaxpyReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_3 )
	{
	    testEvenOddSaxpyReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_1 )
	{
	    testEvenOddSaxpyArgReal( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_2 )
	{
	    testEvenOddSaxpyArgReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_ARG_REAL_EO_3 )
	{
	    testEvenOddSaxpyArgReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_1 )
	{
	    testEvenOddSaxpyVecReal( LatticeExtents{ns4,nt4}, ComplexNumbers {{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_2 )
	{
	    testEvenOddSaxpyVecReal( LatticeExtents{ns8,nt4}, ComplexNumbers {{1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_3 )
	{
	    testEvenOddSaxpyVecReal( LatticeExtents{ns4,nt8}, ComplexNumbers {{-1.,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPBY_EO)

	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_1 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_2 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_3 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_4 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_5 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_6 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_7 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_8 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_CPLX_EO_9 )
	{
	    testEvenOddSaxpbyComplex( LatticeExtents {ns8, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_1 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_2 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_3 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_4 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_5 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_6 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_7 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_8 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_CPLX_EO_9 )
	{
	    testEvenOddSaxpbyArgComplex( LatticeExtents {ns8, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_1 )
	{
	    testEvenOddSaxpbyReal( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_2 )
	{
	    testEvenOddSaxpbyReal( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_3 )
	{
	    testEvenOddSaxpbyReal( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_4 )
	{
	    testEvenOddSaxpbyReal( LatticeExtents {ns4, nt12}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_1 )
	{
	    testEvenOddSaxpbyArgReal( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_2 )
	{
	    testEvenOddSaxpbyArgReal( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_3 )
	{
	    testEvenOddSaxpbyArgReal( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_REAL_EO_4 )
	{
	    testEvenOddSaxpbyArgReal( LatticeExtents {ns4, nt12}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_1 )
	{
	    testEvenOddSaxpbyVecReal( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_2 )
	{
	    testEvenOddSaxpbyVecReal( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_3 )
	{
	    testEvenOddSaxpbyVecReal( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_VEC_4 )
	{
	    testEvenOddSaxpbyVecReal( LatticeExtents {ns4, nt12}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPBYPZ_EO)

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_1 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_2 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_3 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_4 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_5 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_6 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_7 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_8 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_9 )
	{
	    testEvenOddSaxpbypzComplex( LatticeExtents {ns8, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_1 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns4, nt4}, ComplexNumbers {{0.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_2 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns8, nt4}, ComplexNumbers {{1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_3 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns4, nt8}, ComplexNumbers {{0.,1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_4 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns8, nt8}, ComplexNumbers {{0.,-1.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_5 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns12, nt4}, ComplexNumbers {{-1.,0.},{0.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_6 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns4, nt12}, ComplexNumbers {{0.,0.},{-1.,0.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_7 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns12, nt12}, ComplexNumbers {{0.,0.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_8 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns16, nt8}, ComplexNumbers {{0.,1.},{0.,-1.}});
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_9 )
	{
	    testEvenOddSaxpbypzArgComplex( LatticeExtents {ns8, nt16}, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
	testEvenOddGaussianSpinorfield( LatticeExtents{ns12,nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_VEC_AND_SQNORM)

	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_1 )
	{
		testEvenOddSaxVecAndSqnorm( LatticeExtents {ns4,nt4}, ComplexNumbers {{0.,0.}}, 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_2 )
	{
		testEvenOddSaxVecAndSqnorm(  LatticeExtents {ns4,nt4}, ComplexNumbers {{1.,0.}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_3 )
	{
		testEvenOddSaxVecAndSqnorm(  LatticeExtents {ns4,nt4}, ComplexNumbers {{1.,0.1}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_4 )
	{
		testEvenOddSaxVecAndSqnorm(  LatticeExtents {ns4,nt4}, ComplexNumbers {{0.,1.}}, 5);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_5 )
	{
		testEvenOddSaxVecAndSqnorm(  LatticeExtents {ns4,nt4}, ComplexNumbers {{0.01,0.025}}, 15);
	}

BOOST_AUTO_TEST_SUITE_END()
