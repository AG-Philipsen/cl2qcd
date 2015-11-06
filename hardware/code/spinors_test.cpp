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
#include "../../host_functionality/logger.hpp"

#include "mockups.hpp"

static hmc_float sumOfIntegers(const int start, const int end, const int increment) noexcept
{
	// One could also implement some variant of Faulhabers Formula here to save the loop
	hmc_float sum = 0.;
	for(int iteration = start; iteration <= end; iteration += increment)
	{
		sum += iteration;
	}
	return sum;
}

static hmc_float sumOfIntegersSquared(const int end) noexcept
{
	return (2*end*end*end + 3*end*end + end) / 6.; // Faulhaber`s formula
}

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

const ReferenceValues calculateReferenceValues_coldOrZero(const std::string kernelName, const bool isEvenOdd)
{
	if(kernelName == "cold")
	{
		return (isEvenOdd) ? ReferenceValues{0.5} : ReferenceValues{1.};
	}
	if(kernelName == "zero")
	{
		return ReferenceValues{0.};
	}
	else
	{
		return defaultReferenceValues();
	}
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

const ReferenceValues calculateReferenceValues_convert_eo(const int latticeVolume)
{
	return ReferenceValues {latticeVolume * 12., 0.};
}

struct LinearCombinationTestParameters : public SpinorTestParameters
{
	LinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn, const size_t numberOfSpinorsIn, const bool isEvenOddIn):
		SpinorTestParameters(referenceValuesIn, nsIn, ntIn, fillTypesIn, isEvenOddIn), complexNumbers(coefficientsIn), numberOfSpinors(numberOfSpinorsIn){}
	LinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const size_t numberOfSpinorsIn, const bool isEvenOddIn, const int typeOfComparisonIn):
		SpinorTestParameters(referenceValuesIn, nsIn, ntIn, SpinorFillTypes{SpinorFillType::one}, isEvenOddIn, typeOfComparisonIn), complexNumbers(ComplexNumbers{}), numberOfSpinors(numberOfSpinorsIn){}
	LinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const bool isEvenOddIn):
		SpinorTestParameters(referenceValuesIn, nsIn, ntIn, fillTypesIn, isEvenOddIn), complexNumbers(ComplexNumbers {{1.,0.}}), numberOfSpinors(1){}

	const ComplexNumbers complexNumbers;
	const size_t numberOfSpinors;
};

class LinearCombinationTester: public SpinorTester
{
public:
	LinearCombinationTester(const std::string kernelName, const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const LinearCombinationTestParameters testParameters):
			SpinorTester(kernelName, hardwareParameters, kernelParameters, testParameters){

				for (auto coefficient : testParameters.complexNumbers)
				{
					complexNums.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
					complexNums.back()->load(&coefficient);
				}

				for( size_t number = 0; number < testParameters.numberOfSpinors ; number ++)
				{
					if (testParameters.isEvenOdd)
					{
						spinorfieldsEvenOdd.push_back(new hardware::buffers::Spinor(spinorfieldEvenOddElements, device));
						(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfieldsEvenOdd.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfieldsEvenOdd.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
					}
					else
					{
						spinorfields.push_back(new hardware::buffers::Plain<spinor>(spinorfieldElements, device));
						(testParameters.fillTypes.size() < testParameters.numberOfSpinors) ? spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0))) : spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(number)));
					}
				}
			}
protected:
	std::vector<const hardware::buffers::Plain<hmc_complex> *> complexNums;
	std::vector<const hardware::buffers::Plain<spinor> *> spinorfields;
	std::vector<const hardware::buffers::Spinor *> spinorfieldsEvenOdd;
	const hardware::buffers::Plain<spinor> * getOutSpinor() const
	{
		return spinorfields.back();
	}
	const hardware::buffers::Spinor * getOutSpinorEvenOdd() const
	{
		return spinorfieldsEvenOdd.back();
	}
};

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	TesterClass(hardwareParameters, kernelParameters, parametersForThisTest);
}

struct ParameterCollection
{
	ParameterCollection(const hardware::HardwareParametersInterface & hardwareParametersIn, const hardware::code::OpenClKernelParametersInterface & kernelParametersIn):
		hardwareParameters(hardwareParametersIn), kernelParameters(kernelParametersIn) {};
	const hardware::HardwareParametersInterface & hardwareParameters;
	const hardware::code::OpenClKernelParametersInterface & kernelParameters;
};

template<typename TesterClass, typename ParameterClass> void callTest2(const ParameterClass parametersForThisTest)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, fillTypesIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(const int nsIn, const int ntIn, const ComplexNumbers alphaIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, SpinorFillTypes{SpinorFillType::ascendingComplex}, alphaIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(const int nsIn, const int ntIn, const std::string kernelNameIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, kernelNameIn);
	callTest<TesterClass, ParameterClass>(parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest(const int nsIn, const int ntIn, const int typeOfComparisonIn, const bool isEvenOddIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, typeOfComparisonIn, isEvenOddIn);
	callTest2<TesterClass, ParameterClass>(parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE(SPINORTESTER_BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_1") );
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
		BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_2") );
	}

	BOOST_AUTO_TEST_CASE( BUILDFROMPARAMETERS )
	{
		const hardware::HardwareParametersMockup hardwareParameters(4,4);
		const hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(4,4);
		const SpinorTestParameters testParameters;
		BOOST_CHECK_NO_THROW( SpinorTester( "build all kernels", hardwareParameters, kernelParameters, testParameters) );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GLOBAL_SQUARENORM)

	struct SquarenormTestParameters: public LinearCombinationTestParameters
	{
		SquarenormTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn) :
			LinearCombinationTestParameters{ReferenceValues{calculateReferenceValues_globalSquarenorm( calculateSpinorfieldSize(nsIn, ntIn), fillTypesIn)} , nsIn, ntIn, fillTypesIn, false} {};
	};

	struct SquarenormTester: public LinearCombinationTester
	{
		SquarenormTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SquarenormTestParameters & testParameters):
					LinearCombinationTester("global squarenorm", hardwareParameters, kernelParameters, testParameters)
		{
			calcSquarenormAndStoreAsKernelResult(getOutSpinor());
		}
	};

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_1 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>( ns4, nt4, SpinorFillTypes{ SpinorFillType::one} );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_2 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>(ns8, nt4, SpinorFillTypes{ SpinorFillType::ascendingComplex} );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_1 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>(ns8, nt12, SpinorFillTypes{ SpinorFillType::one } );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_2 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>(ns12, nt16, SpinorFillTypes{ SpinorFillType::one } );
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_3 )
	{
		performTest<SquarenormTester, SquarenormTestParameters>(ns16, nt8, SpinorFillTypes{ SpinorFillType::one} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GLOBAL_SQUARENORM_EO)

	struct SquarenormEvenOddTestParameters: public LinearCombinationTestParameters
	{
		SquarenormEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn) :
			LinearCombinationTestParameters{ReferenceValues{calculateReferenceValues_globalSquarenorm( calculateEvenOddSpinorfieldSize(nsIn, ntIn), fillTypesIn)} , nsIn, ntIn, fillTypesIn, true} {};
	};

	struct SquarenormEvenOddTester: public LinearCombinationTester
	{
		SquarenormEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SquarenormEvenOddTestParameters & testParameters):
					LinearCombinationTester("global_squarenorm_eo", hardwareParameters, kernelParameters, testParameters)
		{
			calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
		}
	};
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (ns16, nt8, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (ns16, nt8, SpinorFillTypes{ SpinorFillType::ascendingComplex });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (ns4, nt16, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (ns8, nt4, SpinorFillTypes{ SpinorFillType::one });
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		performTest<SquarenormEvenOddTester, SquarenormEvenOddTestParameters> (ns16, nt16, SpinorFillTypes { SpinorFillType::one });
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	struct ScalarProductTestParameters : public LinearCombinationTestParameters
	{
		ScalarProductTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn):
			LinearCombinationTestParameters(calculateReferenceValues_scalarProduct( calculateSpinorfieldSize(nsIn, ntIn), fillTypesIn), nsIn, ntIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2, false){};
	};

	class ScalarProductTester: public LinearCombinationTester
	{
	public:
		ScalarProductTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ScalarProductTestParameters testParameters):
			LinearCombinationTester("scalar_product", hardwareParameters, kernelParameters, testParameters)
		{
			code->set_complex_to_scalar_product_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0));
			hmc_complex resultTmp;
			complexNums.at(0)->dump(&resultTmp);

			kernelResult[0] = resultTmp.re;
			kernelResult[1] = resultTmp.im;
		}
	};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns4, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns4, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_3 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns4, nt4, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_4 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns4, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns8, nt12, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns8, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		performTest< ScalarProductTester, ScalarProductTestParameters> (ns8, nt16, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	struct ScalarProductEvenOddTestParameters : public LinearCombinationTestParameters
	{
		ScalarProductEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn):
			LinearCombinationTestParameters(calculateReferenceValues_scalarProduct(calculateEvenOddSpinorfieldSize(nsIn, ntIn), fillTypesIn), nsIn, ntIn, fillTypesIn, ComplexNumbers {{1.,0.}}, 2, true) {};
	};

	class ScalarProductEvenOddTester: public LinearCombinationTester
	{
	public:
		ScalarProductEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ScalarProductEvenOddTestParameters testParameters):
			LinearCombinationTester("scalar_product_eo", hardwareParameters, kernelParameters, testParameters)
		{
			code->set_complex_to_scalar_product_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(1), complexNums.at(0));
			hmc_complex resultTmp;
			complexNums.at(0)->dump(&resultTmp);

			kernelResult.at(0) = resultTmp.re;
			kernelResult.at(1) = resultTmp.im;
		}
	};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns8, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns4, nt12, SpinorFillTypes{SpinorFillType::one, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_3 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns8, nt8, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_4 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns4, nt4, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns16, nt4, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns16, nt8, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		performTest<ScalarProductEvenOddTester, ScalarProductEvenOddTestParameters> (ns4, nt16, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO)

	struct ColdAndZeroTestParameters : public LinearCombinationTestParameters
	{
		const std::string kernelName;

		ColdAndZeroTestParameters(const int nsIn, const int ntIn, const std::string kernelNameIn):
			LinearCombinationTestParameters(calculateReferenceValues_coldOrZero(kernelNameIn, false), nsIn, ntIn, SpinorFillTypes{SpinorFillType::one}, false), kernelName(kernelNameIn) {};
	};

	class ColdTester: public LinearCombinationTester
	{
	public:
		ColdTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ColdAndZeroTestParameters testParameters):
			LinearCombinationTester("cold", hardwareParameters, kernelParameters, testParameters)
			{
				code->set_spinorfield_cold_device(getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	class ZeroTester: public LinearCombinationTester
	{
	public:
		ZeroTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ColdAndZeroTestParameters testParameters):
			LinearCombinationTester("zero", hardwareParameters, kernelParameters, testParameters)
			{
				code->set_zero_spinorfield_device(getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		performTest<ColdTester, ColdAndZeroTestParameters> (ns4, nt4, "cold");
	}

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		performTest<ZeroTester, ColdAndZeroTestParameters> (ns4, nt4, "zero");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO_EO)

	struct ColdAndZeroEvenOddTestParameters : public LinearCombinationTestParameters
	{
		const std::string kernelName;

		ColdAndZeroEvenOddTestParameters(const int nsIn, const int ntIn, const std::string kernelNameIn):
			LinearCombinationTestParameters(calculateReferenceValues_coldOrZero(kernelNameIn, true), nsIn, ntIn, SpinorFillTypes{SpinorFillType::one}, true), kernelName(kernelNameIn) {};
	};

	class ColdEvenOddTester: public LinearCombinationTester
	{
	public:
		ColdEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ColdAndZeroEvenOddTestParameters testParameters):
			LinearCombinationTester ("cold_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->set_eoprec_spinorfield_cold_device(getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	class ZeroEvenOddTester: public LinearCombinationTester
	{
	public:
		ZeroEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ColdAndZeroEvenOddTestParameters testParameters):
			LinearCombinationTester("zero_eo",hardwareParameters, kernelParameters, testParameters)
			{
				code->set_zero_spinorfield_eoprec_device(getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};
	
	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		performTest<ColdEvenOddTester, ColdAndZeroEvenOddTestParameters> (ns4, nt4, "cold");
	}

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		performTest<ZeroEvenOddTester, ColdAndZeroEvenOddTestParameters> (ns4, nt4, "zero");
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX)

	struct SaxTestParameters : public LinearCombinationTestParameters
	{
		SaxTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_sax(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 2, false){}
	};

	class SaxTester: public LinearCombinationTester
	{
	public:
		SaxTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxTestParameters testParameters):
			LinearCombinationTester("sax", hardwareParameters, kernelParameters, testParameters)
			{
				code->sax_device(spinorfields.at(0), complexNums.at(0), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		performTest<SaxTester, SaxTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		performTest<SaxTester, SaxTestParameters> (ns8, nt4, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		performTest<SaxTester, SaxTestParameters> (ns4, nt8, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		performTest<SaxTester, SaxTestParameters> (ns16, nt8, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	struct SaxEvenOddTestParameters : public LinearCombinationTestParameters
	{
		SaxEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_sax(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 2, true){}
	};

	class SaxEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxEvenOddTestParameters testParameters):
			LinearCombinationTester("sax_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->sax_eoprec_device(spinorfieldsEvenOdd.at(0), complexNums.at(0), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_EO_1 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_2 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns4, nt8, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_3 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns8, nt4, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_4 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns16, nt8, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	struct SaxpyTestParameters : public LinearCombinationTestParameters
	{
		SaxpyTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxpy(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 3, false){}
	};

	class SaxpyTester: public LinearCombinationTester
	{
	public:
	SaxpyTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyTestParameters testParameters):
			LinearCombinationTester("saxpy", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_device(spinorfields.at(0), spinorfields.at(1), complexNums.at(0), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	class SaxpyArgTester: public LinearCombinationTester
	{
	public:
	SaxpyArgTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyTestParameters testParameters):
			LinearCombinationTester("saxpy_arg", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_device(spinorfields.at(0), spinorfields.at(1), testParameters.complexNumbers.at(0), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns8, nt4, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns4, nt8, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns8, nt8, ComplexNumbers {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns12, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns4, nt12, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns16, nt8, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns8, nt16, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	struct SaxpyEvenOddTestParameters : public LinearCombinationTestParameters
	{
		SaxpyEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxpy(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 3, true){}
	};

	class SaxpyEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxpyEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyEvenOddTestParameters testParameters):
			LinearCombinationTester("saxpy_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(0), complexNums.at(0), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	class SaxpyArgEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxpyArgEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyEvenOddTestParameters testParameters):
			LinearCombinationTester("saxpy_arg_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(0), testParameters.complexNumbers.at(0), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};
	BOOST_AUTO_TEST_CASE( SAXPY_EO_1 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_2 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt4, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_3 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt8, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_4 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt8, ComplexNumbers {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_5 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns12, nt4, ComplexNumbers {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_6 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt12, ComplexNumbers {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_7 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns16, nt8, ComplexNumbers {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_8 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt16, ComplexNumbers {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ)

	struct SaxsbypzTestParameters : public LinearCombinationTestParameters
	{
		SaxsbypzTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxsbypz(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 4, false){}
	};

	class SaxsbypzTester: public LinearCombinationTester
	{
	public:
	SaxsbypzTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzTestParameters testParameters):
			LinearCombinationTester("saxsbypz", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxsbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), complexNums.at(0), complexNums.at(1), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_1 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt4, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt8, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt8, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns12, nt4, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt12, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns12, nt12, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns16, nt8, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt16, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns16, nt16, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)

	struct SaxsbypzEvenOddTestParameters : public LinearCombinationTestParameters
	{
		SaxsbypzEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const ComplexNumbers coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxsbypz(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, 4, true){}
	};

	class SaxsbypzEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxsbypzEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzEvenOddTestParameters testParameters):
				LinearCombinationTester("saxsbypz_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxsbypz_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(1), spinorfieldsEvenOdd.at(2), complexNums.at(0), complexNums.at(1), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt4, ComplexNumbers {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt4, ComplexNumbers {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt8, ComplexNumbers {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt8, ComplexNumbers {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns12, nt4, ComplexNumbers {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt12, ComplexNumbers {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns12, nt12, ComplexNumbers {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns16, nt8, ComplexNumbers {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt16, ComplexNumbers {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(CONVERT_EO)

struct ConvertEvenOddTestParameters: public SpinorTestParameters
{
	const std::string kernelName;

	ConvertEvenOddTestParameters(const int nsIn, const int ntIn, const std::string kernelNameIn):
		SpinorTestParameters(calculateReferenceValues_convert_eo(calculateEvenOddSpinorfieldSize(nsIn, ntIn)), nsIn, ntIn, true),
		kernelName(kernelNameIn) {}
};

	class ConvertEvenOddToOrFromTester: public SpinorTester
	{
	public:
		ConvertEvenOddToOrFromTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ConvertEvenOddTestParameters testParameters):
		SpinorTester("convert_eo_toOrFrom", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);

				in.load( createSpinorfieldWithOnesAndZerosDependingOnSiteParity() );
				code->convert_to_eoprec_device(&in2, &in3, &in) ;

				code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
				doubleBuffer->dump(&kernelResult[1]);
			}
	};

	class ConvertEvenOddTester: public SpinorTester
	{
	public:
		ConvertEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ConvertEvenOddTestParameters testParameters):
			SpinorTester("convert_eo", hardwareParameters, kernelParameters, testParameters)
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

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
		performTest<ConvertEvenOddToOrFromTester, ConvertEvenOddTestParameters> (ns4, nt4, "convert_eo_toOrFrom");
	}

//	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
//	{
//		performTest<ConvertEvenOddToOrFromTester, ConvertEvenOddTestParameters> (ns4, nt4, "convert_eo_toOrFrom");
////		ConvertEvenOddToOrFromTester tester("convert_eo_input_2");
//	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		performTest<ConvertEvenOddTester, ConvertEvenOddTestParameters>(ns4, nt4, "convert_eo");
	}

//	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
//	{
//		ConvertEvenOddTester tester("convert_eo_input_2");
//	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	struct PrngTestParameters : public LinearCombinationTestParameters
	{
	PrngTestParameters(const int nsIn, const int ntIn, const int typeOfComparisonIn, const bool isEvenOddIn) :
		 LinearCombinationTestParameters(calculateReferenceValues_gaussian(), nsIn, ntIn, 1, isEvenOddIn, typeOfComparisonIn), iterations(1000) {};

		const unsigned int iterations;
	};

	class PrngTester: public LinearCombinationTester
	{
	public:
		PrngTester(const std::string kernelName, const ParameterCollection parameterCollection, const PrngTestParameters & testParameters):
					LinearCombinationTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParameters),
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

	class GaussianTester: public PrngTester
	{
	public:
		GaussianTester(const std::string kernelName, const ParameterCollection parameterCollection, const PrngTestParameters & testParameters, const int numberOfElements):
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

		void calculateMean()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				mean += count_sf(&hostOutput[i * numberOfElements], numberOfElements);
			}
			mean /= testParameters.iterations * numberOfElements * 24;
		}

		void calculateVariance()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				variance += calc_var_sf(&hostOutput[i * numberOfElements], numberOfElements, mean);
			}
			variance /= testParameters.iterations * numberOfElements * 24;
		}
	protected:
		const int numberOfElements;
		double mean, variance;
		std::vector<spinor> hostOutput;
		const PrngTestParameters & testParameters;
	};

	class GaussianSpinorfieldTester: public GaussianTester
	{
	public:
		GaussianSpinorfieldTester(const ParameterCollection parameterCollection, const PrngTestParameters testParameters):
					GaussianTester("generate_gaussian_spinorfield", parameterCollection, testParameters, testParameters.getSpinorfieldSize() ){}
		~GaussianSpinorfieldTester()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				code->generate_gaussian_spinorfield_device(this->getOutSpinor(), prngStates);
				this->getOutSpinor()->dump(&hostOutput[i * numberOfElements]);
			}
		}
	};

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_NONEO, 1)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_NONEO )
	{
		performTest<GaussianSpinorfieldTester, PrngTestParameters>(ns8, nt4, comparisonTypes::smallerThan, false);
	}

	class GaussianSpinorfieldEvenOddTester: public GaussianTester
	{
	public:
		GaussianSpinorfieldEvenOddTester(const ParameterCollection parameterCollection, const PrngTestParameters testParameters):
					GaussianTester("generate_gaussian_spinorfield_eo", parameterCollection, testParameters, testParameters.getEvenOddSpinorfieldSize() ) {};
		~GaussianSpinorfieldEvenOddTester()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				code->generate_gaussian_spinorfield_eo_device(this->getOutSpinorEvenOdd(), prngStates);
				this->getOutSpinorEvenOdd()->dump(&hostOutput[i * numberOfElements]);
			}
		}
	};

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GAUSSIAN_EO, 1)

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO )
	{
		performTest<GaussianSpinorfieldEvenOddTester, PrngTestParameters>(ns12, nt4, comparisonTypes::smallerThan, true);
	}

BOOST_AUTO_TEST_SUITE_END()
