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
#include "PrngSpinorTester.hpp"

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
	return ReferenceValues{ 0., 0.707106781};
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

struct LinearCombinationTestParameters : public SpinorTestParameters
{
	LinearCombinationTestParameters(const LatticeExtents lE, const SpinorFillTypes sF) :
		TestParameters(lE), SpinorTestParameters(lE, sF), complexCoefficients(ComplexNumbers{{1.,0.}}), numberOfSpinors(1) {};
	LinearCombinationTestParameters(const LatticeExtents lE, const SpinorFillTypes sF, const ComplexNumbers cN, const size_t numberOfSpinorsIn) :
		TestParameters(lE), SpinorTestParameters(lE, sF), complexCoefficients(cN), numberOfSpinors(numberOfSpinorsIn) {};
	const ComplexNumbers complexCoefficients;
	const NumberOfSpinors numberOfSpinors;
};

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest, const bool needEvenOdd)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, needEvenOdd);
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

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const SpinorFillTypes sF, const ComplexNumbers alphaIn, const int numberOfSpinors, const bool needEvenOdd )
{
	LinearCombinationTestParameters parametersForThisTest(latticeExtendsIn, sF, alphaIn, numberOfSpinors);
	callTest<TesterClass, LinearCombinationTestParameters>(parametersForThisTest, needEvenOdd);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn, const bool fillEvenSitesIn )
{
	SpinorTestParameters parametersForThisTest(latticeExtendsIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest, fillEvenSitesIn);
}

template<typename TesterClass> void performTest(LatticeExtents latticeExtendsIn)
{
	SpinorTestParameters parametersForThisTest(latticeExtendsIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<typename TesterClass, typename ParametersClass>
void performTest(const LatticeExtents latticeExtendsIn, const int iterations, const bool needEvenOdd )
{
	ParametersClass parametersForThisTest(latticeExtendsIn, iterations);
	callTest<TesterClass, PrngSpinorTestParameters>(parametersForThisTest, needEvenOdd);
}

//@todo: this should take the coefficients from another class or be a template class!
struct NonEvenOddLinearCombinationTester : public NonEvenOddSpinorTester
{
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		NonEvenOddSpinorTester(kernelName, pC, tP, rV(calculateSpinorfieldSize(tP.latticeExtents), tP.complexCoefficients))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	NonEvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		NonEvenOddSpinorTester(kernelName, pC, tP, rV)
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
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
		for (auto coefficient : testParameters.complexCoefficients)
		{
			complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
			complexNums.back()->load(&coefficient);
		}
	}
	void loadSpinorfields(const LinearCombinationTestParameters & tP)
	{
		for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
		{
			NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
			spinorfields.push_back(new hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), device));
			(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(sf.createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(sf.createSpinorfield(tP.fillTypes.at(number)));
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

struct EvenOddLinearCombinationTester : public EvenOddSpinorTester
{
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const SpinorFillTypes)):
		EvenOddSpinorTester(kernelName, pC, tP, rV( calculateEvenOddSpinorfieldSize(tP.latticeExtents), tP.fillTypes))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues (*rV) (const int, const ComplexNumbers)):
		EvenOddSpinorTester(kernelName, pC, tP, rV( calculateEvenOddSpinorfieldSize(tP.latticeExtents), tP.complexCoefficients))
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
	}
	EvenOddLinearCombinationTester(const std::string kernelName, const ParameterCollection pC, const LinearCombinationTestParameters tP, const ReferenceValues rV):
		EvenOddSpinorTester(kernelName, pC, tP, rV)
	{
		loadCoefficients(tP);
		loadSpinorfields(tP);
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
			for (auto coefficient : testParameters.complexCoefficients)
			{
				complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
				complexNums.back()->load(&coefficient);
			}
		}
		void loadSpinorfields(const LinearCombinationTestParameters & tP)
		{
			for( size_t number = 0; number < tP.numberOfSpinors ; number ++)
			{
				EvenOddSpinorfieldCreator sf(tP.latticeExtents);
				spinorfields.push_back(new hardware::buffers::Spinor(tP.latticeExtents, device));
				(tP.fillTypes.size() < tP.numberOfSpinors) ? spinorfields.back()->load(sf.createSpinorfield(tP.fillTypes.at(0))) : spinorfields.back()->load(sf.createSpinorfield(tP.fillTypes.at(number)));
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

struct NonEvenGaussianSpinorfieldTester: public PrngSpinorTester
{
	NonEvenGaussianSpinorfieldTester(const ParameterCollection parameterCollection, const PrngSpinorTestParameters testParameters):
		PrngSpinorTester("generate_gaussian_spinorfield", parameterCollection, testParameters, calculateSpinorfieldSize(testParameters.latticeExtents), calculateReferenceValues_gaussian() ){}
	~NonEvenGaussianSpinorfieldTester()
	{
		const hardware::buffers::Plain<spinor> outSpinor(numberOfElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			code->generate_gaussian_spinorfield_device(&outSpinor, prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct EvenOddGaussianSpinorfieldEvenOddTester: public PrngSpinorTester
{
	EvenOddGaussianSpinorfieldEvenOddTester(const ParameterCollection parameterCollection, const PrngSpinorTestParameters testParameters):
				PrngSpinorTester("generate_gaussian_spinorfield_eo", parameterCollection, testParameters, calculateEvenOddSpinorfieldSize(testParameters.latticeExtents), calculateReferenceValues_gaussian() ) {};
	~EvenOddGaussianSpinorfieldEvenOddTester()
	{
		const hardware::buffers::Spinor outSpinor(numberOfElements, device);
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			code->generate_gaussian_spinorfield_eo_device(&outSpinor, prngStates);
			outSpinor.dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct ConvertToEvenOddTester: public SpinorTester
{
	ConvertToEvenOddTester(const ParameterCollection & parameterCollection, const SpinorTestParameters tP, const bool fillEvenSitesIn):
		SpinorTester("convert_to_eo", parameterCollection, tP,
				calculateReferenceValues_convert_eo(calculateEvenOddSpinorfieldSize(tP.latticeExtents), fillEvenSitesIn) )
		{
			const hardware::buffers::Plain<spinor> in(calculateSpinorfieldSize(tP.latticeExtents), device);
			const hardware::buffers::Spinor in2(tP.latticeExtents, device);
			const hardware::buffers::Spinor in3(tP.latticeExtents, device);

			NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
			in.load( sf.createSpinorfieldWithOnesAndZerosDependingOnSiteParity( fillEvenSitesIn ) );
			code->convert_to_eoprec_device(&in2, &in3, &in) ;

			code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);
			code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
			doubleBuffer->dump(&kernelResult[1]);
		}
};

struct ConvertFromEvenOddTester: public SpinorTester
{
	ConvertFromEvenOddTester(const ParameterCollection & parameterCollection, const SpinorTestParameters tP):
		SpinorTester("convert_to_eo", parameterCollection, tP,
				calculateReferenceValues_convertFromEvenOdd(calculateSpinorfieldSize(tP.latticeExtents)) )
		{
			const hardware::buffers::Plain<spinor> in(calculateSpinorfieldSize(tP.latticeExtents), device);
			const hardware::buffers::Spinor in2(tP.latticeExtents, device);
			const hardware::buffers::Spinor in3(tP.latticeExtents, device);

			EvenOddSpinorfieldCreator sf(tP.latticeExtents);
			sf.fillTwoSpinorBuffersDependingOnParity(&in2, &in3);
			code->convert_from_eoprec_device(&in2, &in3, &in);

			code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
			doubleBuffer->dump(&kernelResult[0]);
			}
};

struct SaxsbypzEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxsbypzEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxsbypz_eo", parameterCollection, testParameters, calculateReferenceValues_saxsbypz)
		{
			code->saxsbypz_eoprec_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct SaxsbypzTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxsbypzTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxsbypz", parameterCollection, testParameters, calculateReferenceValues_saxsbypz)
		{
			code->saxsbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
		}
};

struct SquarenormTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SquarenormTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters & testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("global squarenorm", parameterCollection, testParameters, calculateReferenceValues_globalSquarenorm) {}
};

struct SquarenormEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SquarenormEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters & testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("global_squarenorm_eo", parameterCollection, testParameters, calculateReferenceValues_globalSquarenorm) {}
};

struct ScalarProductTester: public NonEvenOddLinearCombinationTester
{
	ScalarProductTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTester("scalar_product", parameterCollection, testParameters, calculateReferenceValues_scalarProduct(calculateSpinorfieldSize(testParameters.latticeExtents), testParameters.fillTypes))
	{
		code->set_complex_to_scalar_product_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
		hmc_complex resultTmp;
		complexNums.at(0)->dump(&resultTmp);

		kernelResult.at(0) = resultTmp.re;
		kernelResult.at(1) = resultTmp.im;
	}
};

struct ScalarProductEvenOddTester: public EvenOddLinearCombinationTester
{
	ScalarProductEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTester("scalar_product_eo", parameterCollection, testParameters, calculateReferenceValues_scalarProduct(calculateEvenOddSpinorfieldSize(testParameters.latticeExtents), testParameters.fillTypes))
	{
		code->set_complex_to_scalar_product_eoprec_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
		hmc_complex resultTmp;
		complexNums.at(0)->dump(&resultTmp);

		kernelResult.at(0) = resultTmp.re;
		kernelResult.at(1) = resultTmp.im;
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
			code->set_spinorfield_cold_device(getOutSpinor());
		}
};

struct ZeroEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ZeroEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult ("zero_eo", parameterCollection, testParameters, calculateReferenceValues_zero())
		{
			code->set_zero_spinorfield_eoprec_device(getOutSpinor());
		}
};

struct ColdEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	ColdEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("cold_eo",parameterCollection, testParameters, calculateReferenceValues_cold(true))
		{
			code->set_eoprec_spinorfield_cold_device(getOutSpinor());
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

struct SaxEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("sax_eo", parameterCollection, testParameters, calculateReferenceValues_sax)
		{
			code->sax_eoprec_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
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

struct SaxpyArgTester: public NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyArgTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		NonEvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_arg", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_device(spinorfields.at(0), spinorfields.at(1), testParameters.complexCoefficients.at(0), getOutSpinor());
		}
};

struct SaxpyEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_eo", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(0), complexNums.at(0), getOutSpinor());
		}
};

struct SaxpyArgEvenOddTester: public EvenOddLinearCombinationTesterWithSquarenormAsKernelResult
{
	SaxpyArgEvenOddTester(const ParameterCollection & parameterCollection, const LinearCombinationTestParameters testParameters):
		EvenOddLinearCombinationTesterWithSquarenormAsKernelResult("saxpy_arg_eo", parameterCollection, testParameters, calculateReferenceValues_saxpy)
		{
			code->saxpy_eoprec_device(spinorfields.at(0), spinorfields.at(0), testParameters.complexCoefficients.at(0), getOutSpinor());
		}
};

void testNonEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<SquarenormTester> (lE, sF, false);
}

void testNonEvenOddScalarProduct( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ScalarProductTester> (lE, sF, ComplexNumbers{{1,0}}, 2, false);
}

void testEvenOddScalarProduct( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ScalarProductEvenOddTester> (lE, sF, ComplexNumbers{{1,0}}, 2, true);
}

void testNonEvenOddSax(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxTester> (lE, cN, 2, false);
}

void testEvenOddSax(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxEvenOddTester> (lE, cN, 2, true);
}

void testNonEvenOddSaxpy(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyTester> (lE, cN, 3, false);
}

void testNonEvenOddSaxpyArg(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyArgTester> (lE, cN, 3, false);
}

void testEvenOddSaxpy(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyEvenOddTester> (lE, cN, 3, true);
}

void testEvenOddSaxpyArg(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxpyArgEvenOddTester> (lE, cN, 3, true);
}

void testNonEvenOddGaussianSpinorfield( const LatticeExtents lE, const int iterations)
{
	performTest<NonEvenGaussianSpinorfieldTester, PrngSpinorTestParameters>(lE, iterations, false);
}

void testEvenOddGaussianSpinorfield( const LatticeExtents lE, const int iterations)
{
	performTest<EvenOddGaussianSpinorfieldEvenOddTester, PrngSpinorTestParameters>(lE, iterations, true);
}

void testConvertToEvenOdd(const LatticeExtents lE, const bool fillEvenSites)
{
	performTest<ConvertToEvenOddTester> (lE, fillEvenSites);
}

void testConvertFromEvenOdd(const LatticeExtents lE)
{
	performTest<ConvertFromEvenOddTester> (lE);
}

void testSaxsbypzEvenOdd( const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxsbypzEvenOddTester> (lE, cN, 4, true);
}

void testNonEvenOddSaxsbypz(const LatticeExtents lE, const ComplexNumbers cN)
{
	performTest<SaxsbypzTester> (lE, cN, 4, false);
}

void testEvenOddSquarenorm( const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<SquarenormEvenOddTester> (lE, sF, true);
}

void testEvenOddCold(const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ColdEvenOddTester> (lE, sF, true);
}

void testEvenOddZero(const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ZeroEvenOddTester> (lE, sF, true);
}

void testNonEvenOddCold(const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ColdTester> (lE, sF, false);
}

void testNonEvenOddZero(const LatticeExtents lE, const SpinorFillTypes sF)
{
	performTest<ZeroTester> (lE, sF, false);
}

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
		testNonEvenOddZero(LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD)

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		testNonEvenOddCold(LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ZERO_EO)

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		testEvenOddZero( LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_EO)

	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		testEvenOddCold(LatticeExtents{ns4, nt4}, SpinorFillTypes{SpinorFillType::one});
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
		testConvertFromEvenOdd(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_NONEO, 2)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_NONEO )
	{
		testNonEvenOddGaussianSpinorfield( LatticeExtents{ns8, nt4}, 1000 );
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_EO, 2)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO )
	{
		testEvenOddGaussianSpinorfield( LatticeExtents{ns12, nt4}, 1000 );
	}

BOOST_AUTO_TEST_SUITE_END()
