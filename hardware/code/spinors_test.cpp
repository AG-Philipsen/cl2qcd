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

static int calculateLatticeVolume(const int nsIn, const int ntIn) noexcept
{
	return 	nsIn * nsIn * nsIn * ntIn;
}

static int calculateSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateLatticeVolume(nsIn, ntIn);
}

static int calculateEvenOddSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateSpinorfieldSize(nsIn, ntIn) / 2;
}

static ReferenceValues defaultReferenceValues()
{
	return ReferenceValues{-1.23456};
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

const ReferenceValues calculateReferenceValues_sax (const int latticeVolume, Coefficients alphaIn)
{
	return ReferenceValues{ (alphaIn.at(0).im * alphaIn.at(0).im + alphaIn.at(0).re * alphaIn.at(0).re) * latticeVolume * sumOfIntegersSquared(24)};
}

const ReferenceValues calculateReferenceValues_saxpy(const int latticeVolume, Coefficients alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, Coefficients {{1. - alphaIn.at(0).re, 0. - alphaIn.at(0).im}}).at(0)};
}

const ReferenceValues calculateReferenceValues_saxsbypz(const int latticeVolume, Coefficients alphaIn)
{
	return ReferenceValues{calculateReferenceValues_sax(latticeVolume, Coefficients {{1. + alphaIn.at(0).re + alphaIn.at(1).re, 0. + alphaIn.at(0).im + alphaIn.at(1).im}}).at(0)};
}

struct LinearCombinationTestParameters : public SpinorTestParameters
{
	Coefficients coefficients;

	LinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn, const bool isEvenOddIn):
		SpinorTestParameters(referenceValuesIn, nsIn, ntIn, fillTypesIn, isEvenOddIn), coefficients(coefficientsIn){}
	LinearCombinationTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const bool isEvenOddIn):
			SpinorTestParameters(referenceValuesIn, nsIn, ntIn, fillTypesIn, isEvenOddIn), coefficients(Coefficients {{1.,0.}}){}
};

class LinearCombinationTester: public SpinorTester
{
public:
	LinearCombinationTester(std::string kernelName, const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, LinearCombinationTestParameters testParameters):
			SpinorTester(kernelName, hardwareParameters, kernelParameters, testParameters){

				for (auto coefficient : testParameters.coefficients)
				{
					coeff.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
					coeff.back()->load(&coefficient);
				}

				for( auto number = 0; number < 4; number ++)
				{
					if (testParameters.isEvenOdd)
					{
						spinorfieldsEvenOdd.push_back(new hardware::buffers::Spinor(spinorfieldEvenOddElements, device));
						spinorfieldsEvenOdd.back()->load(createSpinorfield(testParameters.fillTypes.at(0)));
					}
					else
					{
						spinorfields.push_back(new hardware::buffers::Plain<spinor>(spinorfieldElements, device));
						spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0)));
					}
				}
			}
protected:
	std::vector<const hardware::buffers::Plain<hmc_complex> *> coeff;
	std::vector<const hardware::buffers::Plain<spinor> *> spinorfields;
	std::vector<const hardware::buffers::Spinor *> spinorfieldsEvenOdd;
	const hardware::buffers::Plain<spinor> * getOutSpinor()
	{
		return spinorfields.back();
	}
	const hardware::buffers::Spinor * getOutSpinorEvenOdd()
	{
		return spinorfieldsEvenOdd.back();
	}
};

template<typename TesterClass, typename ParameterClass> void performTest( const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, fillTypesIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	TesterClass(hardwareParameters, kernelParameters, parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest( const int nsIn, const int ntIn, const Coefficients alphaIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, SpinorFillTypes{SpinorFillType::ascendingComplex}, alphaIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	TesterClass(hardwareParameters, kernelParameters, parametersForThisTest);
}

template<typename TesterClass, typename ParameterClass> void performTest( const int nsIn, const int ntIn, const std::string kernelNameIn )
{
	ParameterClass parametersForThisTest(nsIn, ntIn, kernelNameIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.isEvenOdd);
	TesterClass(hardwareParameters, kernelParameters, parametersForThisTest);
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
		hardware::HardwareParametersMockup hardwareParameters(4,4);
		hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(4,4);
		SpinorTestParameters testParameters{std::vector<double> {-1.234}, 4,4, SpinorFillTypes{SpinorFillType::one}, false};
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

	struct ScalarProductTestParameters : public SpinorTestParameters
	{
		ScalarProductTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn):
			SpinorTestParameters(calculateReferenceValues_scalarProduct( calculateSpinorfieldSize(nsIn, ntIn), fillTypesIn), nsIn, ntIn, fillTypesIn, false){};
	};

	class ScalarProductTester: public SpinorTester
	{
	public:
		ScalarProductTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ScalarProductTestParameters testParameters):
			SpinorTester("scalar_product", hardwareParameters, kernelParameters, testParameters)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
			in.load(createSpinorfield(testParameters.fillTypes[0]));
			in2.load(createSpinorfield(testParameters.fillTypes[1]));
			hardware::buffers::Plain<hmc_complex> sqnorm(1, device);

			code->set_complex_to_scalar_product_device(&in, &in2, &sqnorm);
			hmc_complex resultTmp;
			sqnorm.dump(&resultTmp);

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

	struct ScalarProductEvenOddTestParameters : public SpinorTestParameters
	{
		ScalarProductEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn):
			SpinorTestParameters(calculateReferenceValues_scalarProduct(calculateEvenOddSpinorfieldSize(nsIn, ntIn), fillTypesIn), nsIn, ntIn, fillTypesIn, true) {};
	};

	class ScalarProductEvenOddTester: public SpinorTester
	{
	public:
		ScalarProductEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const ScalarProductEvenOddTestParameters testParameters):
			SpinorTester("scalar_product_eo", hardwareParameters, kernelParameters, testParameters)
		{
			const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
			const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield (testParameters.fillTypes.at(0)));
			in2.load(createSpinorfield(testParameters.fillTypes.at(1)));
			hardware::buffers::Plain<hmc_complex> sqnorm(1, device);

			code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
			hmc_complex resultTmp;
			sqnorm.dump(&resultTmp);

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
		SaxTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_sax(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, false){}
	};

	class SaxTester: public LinearCombinationTester
	{
	public:
		SaxTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxTestParameters testParameters):
			LinearCombinationTester("sax", hardwareParameters, kernelParameters, testParameters)
			{
				code->sax_device(spinorfields.at(0), coeff.at(0), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		performTest<SaxTester, SaxTestParameters> (ns4, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		performTest<SaxTester, SaxTestParameters> (ns8, nt4, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		performTest<SaxTester, SaxTestParameters> (ns4, nt8, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		performTest<SaxTester, SaxTestParameters> (ns16, nt8, Coefficients {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	struct SaxEvenOddTestParameters : public LinearCombinationTestParameters
	{
		SaxEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_sax(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, true){}
	};

	class SaxEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxEvenOddTestParameters testParameters):
			LinearCombinationTester("sax_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->sax_eoprec_device(spinorfieldsEvenOdd.at(0), coeff.at(0), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_EO_1 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns4, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_2 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns4, nt8, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_3 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns8, nt4, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_4 )
	{
		performTest<SaxEvenOddTester, SaxEvenOddTestParameters> (ns16, nt8, Coefficients {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	struct SaxpyTestParameters : public LinearCombinationTestParameters
	{
		SaxpyTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxpy(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, false){}
	};

	class SaxpyTester: public LinearCombinationTester
	{
	public:
	SaxpyTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyTestParameters testParameters):
			LinearCombinationTester("saxpy", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_device(spinorfields.at(0), spinorfields.at(1), coeff.at(0), getOutSpinor());
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
				code->saxpy_device(spinorfields.at(0), spinorfields.at(1), coeff.at(0), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns4, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns8, nt4, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		performTest<SaxpyTester, SaxpyTestParameters> (ns4, nt8, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns8, nt8, Coefficients {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns12, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns4, nt12, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns16, nt8, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
		performTest<SaxpyArgTester, SaxpyTestParameters> (ns8, nt16, Coefficients {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

struct SaxpyEvenOddTestParameters : public LinearCombinationTestParameters
{
	SaxpyEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
		LinearCombinationTestParameters(calculateReferenceValues_saxpy(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, true){}
};

	class SaxpyEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxpyEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyEvenOddTestParameters testParameters):
			LinearCombinationTester("saxpy_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxpy_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(0), coeff.at(0), getOutSpinorEvenOdd());
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
				code->saxpy_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(0), testParameters.coefficients.at(0), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};
	BOOST_AUTO_TEST_CASE( SAXPY_EO_1 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_2 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt4, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_3 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt8, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_4 )
	{
		performTest<SaxpyEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt8, Coefficients {{1.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_5 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns12, nt4, Coefficients {{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_6 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns4, nt12, Coefficients {{1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_7 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns16, nt8, Coefficients {{0.,1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_8 )
	{
		performTest<SaxpyArgEvenOddTester, SaxpyEvenOddTestParameters> (ns8, nt16, Coefficients {{1.,1.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ)

	struct SaxsbypzTestParameters : public LinearCombinationTestParameters
	{
		SaxsbypzTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
			LinearCombinationTestParameters(calculateReferenceValues_saxsbypz(calculateSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, false){}
	};

	class SaxsbypzTester: public LinearCombinationTester
	{
	public:
	SaxsbypzTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzTestParameters testParameters):
			LinearCombinationTester("saxsbypz", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxsbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), coeff.at(0), coeff.at(1), getOutSpinor());
				calcSquarenormAndStoreAsKernelResult(getOutSpinor());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_1 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt4, Coefficients {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt4, Coefficients {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt8, Coefficients {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt8, Coefficients {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns12, nt4, Coefficients {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns4, nt12, Coefficients {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns12, nt12, Coefficients {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns16, nt8, Coefficients {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns8, nt16, Coefficients {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		performTest<SaxsbypzTester, SaxsbypzTestParameters> (ns16, nt16, Coefficients {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)

struct SaxsbypzEvenOddTestParameters : public LinearCombinationTestParameters
{
	SaxsbypzEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients coefficientsIn):
		LinearCombinationTestParameters(calculateReferenceValues_saxsbypz(calculateEvenOddSpinorfieldSize(nsIn, ntIn), coefficientsIn), nsIn, ntIn, fillTypesIn, coefficientsIn, true){}
};

	class SaxsbypzEvenOddTester: public LinearCombinationTester
	{
	public:
		SaxsbypzEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzEvenOddTestParameters testParameters):
				LinearCombinationTester("saxsbypz_eo", hardwareParameters, kernelParameters, testParameters)
			{
				code->saxsbypz_eoprec_device(spinorfieldsEvenOdd.at(0), spinorfieldsEvenOdd.at(1), spinorfieldsEvenOdd.at(2), coeff.at(0), coeff.at(1), getOutSpinorEvenOdd());
				calcSquarenormEvenOddAndStoreAsKernelResult(getOutSpinorEvenOdd());
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt4, Coefficients {{0.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt4, Coefficients {{1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt8, Coefficients {{0.,1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt8, Coefficients {{0.,-1.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns12, nt4, Coefficients {{-1.,0.},{0.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns4, nt12, Coefficients {{0.,0.},{-1.,0.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns12, nt12, Coefficients {{0.,0.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns16, nt8, Coefficients {{0.,1.},{0.,-1.}});
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
	{
		performTest<SaxsbypzEvenOddTester, SaxsbypzEvenOddTestParameters> (ns8, nt16, Coefficients {{-0.5,0.},{-0.5,0.}});
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	class ConvertEvenOddTester: public SpinorTester
	{
	public:
		ConvertEvenOddTester(std::string inputfile, bool toOrFrom = true):
			SpinorTester("convert_eo", inputfile, 2)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);

				

				if ( toOrFrom )
				{
					in.load( createSpinorfieldWithOnesAndZerosDependingOnSiteParity() );
					code->convert_to_eoprec_device(&in2, &in3, &in) ;

					code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
					doubleBuffer->dump(&kernelResult[0]);
					code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
					doubleBuffer->dump(&kernelResult[1]);
					
					if (evenOrOdd)
					{
						referenceValue[0] = spinorfieldEvenOddElements*12.; 
						referenceValue[1] = 0.;
					}
					else
					{ 
						referenceValue[0] = 0.;
						referenceValue[1] = spinorfieldEvenOddElements*12.; 
					}
				} else
				{
					fillTwoSpinorBuffersDependingOnParity(&in2, &in3);
					code->convert_from_eoprec_device(&in2, &in3, &in);
					
					code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
					doubleBuffer->dump(&kernelResult[0]);
					
					referenceValue[0] = spinorfieldEvenOddElements*12.; 
				}

				

			}
	};

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
		ConvertEvenOddTester tester("convert_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
		ConvertEvenOddTester tester("convert_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		ConvertEvenOddTester tester("convert_eo_input_1", false);
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
		ConvertEvenOddTester tester("convert_eo_input_2", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	class GaussianTester: public SpinorTester
	{
	public:
		GaussianTester(std::string inputfile):
		  SpinorTester("generate_gaussian_spinorfield", inputfile, 1, 2)
			{
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				spinor * outHost;
				outHost = new spinor[spinorfieldElements * iterations];
				BOOST_REQUIRE(out);
				
				auto prng_buf = prng->get_buffers().at(0);
		
				double sum = 0;
				for (int i = 0; i < iterations; i++) {
					code->generate_gaussian_spinorfield_device(&out, prng_buf);
					out.dump(&outHost[i * spinorfieldElements]);
					sum += count_sf(&outHost[i * spinorfieldElements], spinorfieldElements);
				}
				sum /= iterations * spinorfieldElements * 24;
				kernelResult[0] = sum;

				if(calcVariance) 
				{
					double var = 0.;
					for (int i = 0; i < iterations; i++) {
						var += calc_var_sf(&outHost[i * spinorfieldElements], spinorfieldElements, sum);
					}
					var /= iterations * spinorfieldElements * 24;
					kernelResult[0] =  sqrt(var);
				}
			}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
		GaussianTester tester("gaussian_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_2 )
	{
		GaussianTester tester("gaussian_input_2");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_3 )
	{
		GaussianTester tester("gaussian_input_3");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_4 )
	{
		GaussianTester tester("gaussian_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	class GaussianEvenOddTester: public SpinorTester
	{
	public:
		GaussianEvenOddTester(std::string inputfile):
		  SpinorTester("generate_gaussian_spinorfield_eo", inputfile, 1, 2)
			{
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				spinor * outHost;
				outHost = new spinor[spinorfieldEvenOddElements * iterations];
				BOOST_REQUIRE(out);
				
				auto prng_buf = prng->get_buffers().at(0);

				double sum = 0;
				for (int i = 0; i < iterations; i++) {
					code->generate_gaussian_spinorfield_eo_device(&out, prng_buf);
					out.dump(&outHost[i * spinorfieldEvenOddElements]);
					sum += count_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements);
				}
				sum /= iterations * spinorfieldElements * 24;
				kernelResult[0] = sum;

				if(calcVariance) 
				{
					double var = 0.;
					for (int i = 0; i < iterations; i++) {
						var += calc_var_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements, sum);
					}
					var /= iterations * spinorfieldEvenOddElements * 24;
					kernelResult[0] =  sqrt(var);
				}
			}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_2 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_3 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_4 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()
