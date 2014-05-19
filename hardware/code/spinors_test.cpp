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

BOOST_AUTO_TEST_SUITE(SPINORTESTER_BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_1") );
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
		BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_2") );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GLOBAL_SQUARENORM)

	class SquarenormTester: public SpinorTester
	{
	public:
		SquarenormTester(std::string inputfile):
			SpinorTester("global_squarenorm", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				in.load(createSpinorfield(spinorfieldElements));
				calcSquarenormAndStoreAsKernelResult(&in);
			}
	};

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_1 )
	{
		SquarenormTester("/global_squarenorm_input_1");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_2 )
	{
		SquarenormTester("/global_squarenorm_input_2");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_1 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_2 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_3 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GLOBAL_SQUARENORM_EO)

	class SquarenormEvenOddTester: public SpinorTester
	{
	public:
		SquarenormEvenOddTester(std::string inputfile):
			SpinorTester("global_squarenorm_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				in.load(createSpinorfield(spinorfieldEvenOddElements));
				calcSquarenormEvenOddAndStoreAsKernelResult(&in);
			}
	};
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	class ScalarProductTester: public SpinorTester
	{
	public:
		ScalarProductTester(std::string inputfile):
			SpinorTester("scalar_product", inputfile, 2)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
				in.load(createSpinorfield(spinorfieldElements, 123));
				in2.load(createSpinorfield(spinorfieldElements, 456));
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
		ScalarProductTester scalarProductTester("scalar_product_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		ScalarProductTester scalarProductTester("scalar_product_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	class ScalarProductEvenOddTester: public SpinorTester
	{
	public:
		ScalarProductEvenOddTester(std::string inputfile):
			SpinorTester("scalar_product_eo", inputfile, 2)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
				in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
				hardware::buffers::Plain<hmc_complex> sqnorm(1, device);

				code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
				hmc_complex resultTmp;
				sqnorm.dump(&resultTmp);
				
				kernelResult[0] = resultTmp.re;
				kernelResult[1] = resultTmp.im;
			}
	};
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_3");
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

	class SaxTester: public SpinorTester
	{
	public:
		SaxTester(std::string inputfile):
			SpinorTester("sax", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(spinorfieldElements, 123));
				alpha.load(&alpha_host);

				code->sax_device(&in, &alpha, &out);
				calcSquarenormAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		SaxTester tester ("sax_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		SaxTester tester ("sax_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		SaxTester tester ("sax_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		SaxTester tester ("sax_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAX_5 )
	{
		SaxTester tester ("sax_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAX_6 )
	{
		SaxTester tester ("sax_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAX_7 )
	{
		SaxTester tester ("sax_input_7");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	class SaxEvenOddTester: public SpinorTester
	{
	public:
		SaxEvenOddTester(std::string inputfile):
			SpinorTester("sax_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
				alpha.load(&alpha_host);

				code->sax_eoprec_device(&in, &alpha, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_EO_1 )
	{
		SaxEvenOddTester tester("sax_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_2 )
	{
		SaxEvenOddTester tester("sax_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_3 )
	{
		SaxEvenOddTester tester("sax_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_4 )
	{
		SaxEvenOddTester tester("sax_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_5 )
	{
		SaxEvenOddTester tester("sax_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_6 )
	{
		SaxEvenOddTester tester("sax_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_7 )
	{
		SaxEvenOddTester tester("sax_eo_input_7");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	class SaxpyTester: public SpinorTester
	{
	public:
		SaxpyTester(std::string inputfile, bool switcher):
			SpinorTester("saxpy", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(spinorfieldElements, 123));
				in2.load(createSpinorfield(spinorfieldElements, 456));
				alpha.load(&alpha_host);

				(switcher) ? code->saxpy_device(&in, &in2, &alpha, &out) : code->saxpy_device(&in, &in2, alpha_host, &out);
				calcSquarenormAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		SaxpyTester tester("saxpy_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		SaxpyTester tester("saxpy_input_2", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		SaxpyTester tester("saxpy_input_3", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		SaxpyTester tester("saxpy_input_4", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
		SaxpyTester tester("/saxpy_input_5", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
		SaxpyTester tester("saxpy_input_6", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
		SaxpyTester tester("saxpy_input_7", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
		SaxpyTester tester("saxpy_input_8", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_9 )
	{
		SaxpyTester tester("saxpy_input_9", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_10 )
	{
		SaxpyTester tester("saxpy_input_10", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_11 )
	{
		SaxpyTester tester("saxpy_input_11", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_12 )
	{
		SaxpyTester tester("saxpy_input_12", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_13 )
	{
		SaxpyTester tester("saxpy_input_13", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_14 )
	{
		SaxpyTester tester("saxpy_input_14", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_1 )
	{
		SaxpyTester tester("saxpy_input_1", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_2 )
	{
		SaxpyTester tester("saxpy_input_2", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_3 )
	{
		SaxpyTester tester("saxpy_input_3", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_4 )
	{
		SaxpyTester tester("saxpy_input_4", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_5 )
	{
		SaxpyTester tester("saxpy_input_5", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_6 )
	{
		SaxpyTester tester("saxpy_input_6", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_7 )
	{
		SaxpyTester tester("saxpy_input_7", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_8 )
	{
		SaxpyTester tester("saxpy_input_8", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_9 )
	{
		SaxpyTester tester("saxpy_input_9", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_10 )
	{
		SaxpyTester tester("saxpy_input_10", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_11 )
	{
		SaxpyTester tester("saxpy_input_11", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_12 )
	{
		SaxpyTester tester("saxpy_input_12", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_13 )
	{
		SaxpyTester tester("saxpy_input_13", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_14 )
	{
		SaxpyTester tester("saxpy_input_14", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	class SaxpyEvenOddTester: public SpinorTester
	{
	public:
		SaxpyEvenOddTester(std::string inputfile, bool switcher):
			SpinorTester("saxpy_eo or saxpy_arg_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
				in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
				alpha.load(&alpha_host);

				(switcher) ? code->saxpy_eoprec_device(&in, &in2, &alpha, &out) : code->saxpy_eoprec_device(&in, &in2, alpha_host, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_EO_1 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_2 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_2", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_3 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_3", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_4 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_4", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_5 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_5", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_6 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_6", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_7 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_7", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_8 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_8", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_9 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_9", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_10 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_10", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_11 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_11", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_12 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_12", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_13 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_13", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_14 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_14", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_1 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_1", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_2 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_2", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_3 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_3", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_4 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_4", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_5 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_5", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_6 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_6", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_7 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_7", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_8 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_8", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_9 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_9", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_10 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_10", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_11 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_11", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_12 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_12", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_13 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_13", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_14 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_14", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ)

	class SaxsbypzTester: public SpinorTester
	{
	public:
		SaxsbypzTester(std::string inputfile):
			SpinorTester("saxsbypz", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in2(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> in3(spinorfieldElements, device);
				const hardware::buffers::Plain<spinor> out(spinorfieldElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);

				in.load(createSpinorfield(spinorfieldElements, 123));
				in2.load(createSpinorfield(spinorfieldElements, 456));
				in3.load(createSpinorfield(spinorfieldElements, 789));
				alpha.load(&alpha_host);
				beta.load(&beta_host);

				code->saxsbypz_device(&in, &in2, &in3, &alpha, &beta, &out);
				calcSquarenormAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_1 )
	{
		SaxsbypzTester tester("saxsbypz_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		SaxsbypzTester tester("saxsbypz_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		SaxsbypzTester tester("saxsbypz_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		SaxsbypzTester tester("saxsbypz_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		SaxsbypzTester tester("saxsbypz_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		SaxsbypzTester tester("saxsbypz_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		SaxsbypzTester tester("saxsbypz_input_7");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		SaxsbypzTester tester("saxsbypz_input_8");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		SaxsbypzTester tester("saxsbypz_input_9");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		SaxsbypzTester tester("saxsbypz_input_10");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)

	class SaxsbypzEvenOddTester: public SpinorTester
	{
	public:
		SaxsbypzEvenOddTester(std::string inputfile):
			SpinorTester("saxsbypz_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor in3(spinorfieldEvenOddElements, device);
				const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);

				in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
				in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
				in3.load(createSpinorfield(spinorfieldEvenOddElements, 789));
				alpha.load(&alpha_host);
				beta.load(&beta_host);

				code->saxsbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);
				calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_7");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_8");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_9");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_10 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_10");
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
					var /= iterations * spinorfieldEvenOddElements * 24;
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
