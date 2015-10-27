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

	struct SquarenormTestParameters: public SpinorTestParameters
	{
		SquarenormTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn) :
			SpinorTestParameters{ReferenceValues{calculateReferenceValues_globalSquarenorm( calculateSpinorfieldSize(nsIn, ntIn), fillTypesIn)} , nsIn, ntIn, fillTypesIn, false} {};
	};

	struct SquarenormTester: public SpinorTester
	{
		SquarenormTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SquarenormTestParameters & testParameters):
					SpinorTester("global squarenorm", hardwareParameters, kernelParameters, testParameters)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			in.load(createSpinorfield( testParameters.fillTypes.at(0)) );
			calcSquarenormAndStoreAsKernelResult(&in);
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

	struct SquarenormEvenOddTestParameters: public SpinorTestParameters
	{
		SquarenormEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn) :
			SpinorTestParameters{ReferenceValues{calculateReferenceValues_globalSquarenorm( calculateEvenOddSpinorfieldSize(nsIn, ntIn), fillTypesIn)} , nsIn, ntIn, fillTypesIn, true} {};
	};

	struct SquarenormEvenOddTester: public SpinorTester
	{
		SquarenormEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SquarenormEvenOddTestParameters & testParameters):
					SpinorTester("global_squarenorm_eo", hardwareParameters, kernelParameters, testParameters)
		{
			const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
			in.load( createSpinorfield( testParameters.fillTypes.at(0) ) );
			calcSquarenormEvenOddAndStoreAsKernelResult(&in);
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

	class ColdAndZeroTester: public SpinorTester
	{
	public:
		ColdAndZeroTester(std::string inputfile, bool switcher):
			SpinorTester("cold or zero", inputfile)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				in.load(createSpinorfield(spinorfieldElements));
				(switcher) ? code->set_spinorfield_cold_device(&in) : 	code->set_zero_spinorfield_device(&in);
				calcSquarenormAndStoreAsKernelResult(&in);
			}
	};

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		ColdAndZeroTester tester("cold_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		ColdAndZeroTester tester("zero_input_1", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO_EO)

	class ColdAndZeroEvenOddTester: public SpinorTester
	{
	public:
		ColdAndZeroEvenOddTester(std::string inputfile, bool switcher):
			SpinorTester("cold or zero eo", inputfile)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				in.load(createSpinorfield(spinorfieldEvenOddElements));
				(switcher) ? code->set_eoprec_spinorfield_cold_device(&in) : 	code->set_zero_spinorfield_eoprec_device(&in);
				calcSquarenormEvenOddAndStoreAsKernelResult(&in);
			}
	};
	
	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		ColdAndZeroEvenOddTester tester("cold_eo_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		ColdAndZeroEvenOddTester tester("zero_eo_input_1",  false);
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX)

	struct SaxTestParameters : public SpinorTestParameters
	{
		Coefficients alpha;

		SaxTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
			SpinorTestParameters(calculateReferenceValues_sax(calculateSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, false), alpha(alphaIn){}
	};

	class SaxTester: public SpinorTester
	{
	public:
		SaxTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxTestParameters testParameters):
			SpinorTester("sax", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->sax_device(&in, &alpha, &out);
				calcSquarenormAndStoreAsKernelResult(&out);
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

	struct SaxEvenOddTestParameters : public SpinorTestParameters
	{
		Coefficients alpha;

		SaxEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
			SpinorTestParameters(calculateReferenceValues_sax(calculateEvenOddSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, true), alpha(alphaIn){}
	};

	class SaxEvenOddTester: public SpinorTester
	{
	public:
		SaxEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxEvenOddTestParameters testParameters):
		SpinorTester("sax_eo", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->sax_eoprec_device(&in, &alpha, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
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

	struct SaxpyTestParameters : public SpinorTestParameters
	{
		Coefficients alpha;

		SaxpyTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
			SpinorTestParameters(calculateReferenceValues_saxpy(calculateSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, false), alpha(alphaIn){}
	};

	class SaxpyTester: public SpinorTester
	{
	public:
	SaxpyTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyTestParameters testParameters):
			SpinorTester("saxpy", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in2.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->saxpy_device(&in, &in2, &alpha, &out);
				calcSquarenormAndStoreAsKernelResult(&out);
			}
	};

	class SaxpyArgTester: public SpinorTester
	{
	public:
	SaxpyArgTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyTestParameters testParameters):
			SpinorTester("saxpy_arg", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in2.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->saxpy_device(&in, &in2, testParameters.alpha.at(0), &out);
				calcSquarenormAndStoreAsKernelResult(&out);
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

struct SaxpyEvenOddTestParameters : public SpinorTestParameters
{
	Coefficients alpha;

	SaxpyEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
		SpinorTestParameters(calculateReferenceValues_saxpy(calculateEvenOddSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, true), alpha(alphaIn){}
};

	class SaxpyEvenOddTester: public SpinorTester
	{
	public:
		SaxpyEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyEvenOddTestParameters testParameters):
			SpinorTester("saxpy_eo", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in2.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->saxpy_eoprec_device(&in, &in2, &alpha, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			}
	};

	class SaxpyArgEvenOddTester: public SpinorTester
	{
	public:
		SaxpyArgEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxpyEvenOddTestParameters testParameters):
			SpinorTester("saxpy_arg_eo", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in2.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));

				code->saxpy_eoprec_device(&in, &in2, testParameters.alpha.at(0), &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
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

struct SaxsbypzTestParameters : public SpinorTestParameters
{
	Coefficients alpha;

	SaxsbypzTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
		SpinorTestParameters(calculateReferenceValues_saxsbypz(calculateSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, false), alpha(alphaIn){}
};

	class SaxsbypzTester: public SpinorTester
	{
	public:
	SaxsbypzTester(const hardware::HardwareParametersInterface & hardwareParameters,
			const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzTestParameters testParameters):
			SpinorTester("saxsbypz", hardwareParameters, kernelParameters, testParameters)
			{
				std::vector<const hardware::buffers::Plain<hmc_complex> *> coeff;
				for (auto coefficient : testParameters.alpha)
				{
					coeff.push_back(new hardware::buffers::Plain<hmc_complex>(1, device));
					coeff.back()->load(&coefficient);
				}

				std::vector<const hardware::buffers::Plain<spinor> *> spinorfields;
				for( auto number = 0; number < 4; number ++)
				{
					spinorfields.push_back(new hardware::buffers::Plain<spinor>(spinorfieldElements, device));
					spinorfields.back()->load(createSpinorfield(testParameters.fillTypes.at(0)));
				}

				code->saxsbypz_device(spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), coeff.at(0), coeff.at(1), spinorfields.at(3));
				calcSquarenormAndStoreAsKernelResult(spinorfields.at(3));
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

struct SaxsbypzEvenOddTestParameters : public SpinorTestParameters
{
	Coefficients alpha;

	SaxsbypzEvenOddTestParameters(const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, Coefficients alphaIn):
		SpinorTestParameters(calculateReferenceValues_saxsbypz(calculateEvenOddSpinorfieldSize(nsIn, ntIn), alphaIn), nsIn, ntIn, fillTypesIn, true), alpha(alphaIn){}
};

	class SaxsbypzEvenOddTester: public SpinorTester
	{
	public:
		SaxsbypzEvenOddTester(const hardware::HardwareParametersInterface & hardwareParameters,
				const hardware::code::OpenClKernelParametersInterface & kernelParameters, const SaxsbypzEvenOddTestParameters testParameters):
			SpinorTester("saxsbypz_eo", hardwareParameters, kernelParameters, testParameters)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);

				in.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in2.load(createSpinorfield(testParameters.fillTypes.at(0)));
				in3.load(createSpinorfield(testParameters.fillTypes.at(0)));
				alpha.load(&testParameters.alpha.at(0));
				beta.load(&testParameters.alpha.at(1));

				code->saxsbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
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
