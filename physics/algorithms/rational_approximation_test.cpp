/** @file
 * Tests of the Rational_Approximation algorithm
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "rational_approximation.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::rational_approximation
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>


BOOST_AUTO_TEST_CASE(initialization)
{
	if(boost::unit_test::framework::master_test_suite().argc !=2)
	  logger.error() << "In the initialization test the sourcefile for the approx. is needed!";
	BOOST_REQUIRE_EQUAL(boost::unit_test::framework::master_test_suite().argc, 2);
	
        using namespace physics::algorithms;
	
	Rational_Approximation approx(6,1,2,1e-5,1);
	logger.info() << approx;
	int ord = approx.Get_order();
	hmc_float a0 = approx.Get_a0();
	std::vector<hmc_float> a = approx.Get_a();
	std::vector<hmc_float> b = approx.Get_b();
	
	Rational_Coefficients coeff_alloc(ord);
	Rational_Coefficients coeff(ord, a0, a, b);
	int ord2 = coeff.Get_order();
	hmc_float a02 = coeff.Get_a0();
	std::vector<hmc_float> a2 = coeff.Get_a();
	std::vector<hmc_float> b2 = coeff.Get_b();
	BOOST_CHECK_CLOSE(a0, a02, 1.e-8);
	BOOST_REQUIRE_EQUAL(ord, ord2);
	for(int i=0; i<ord; i++){
		BOOST_CHECK_CLOSE(a[i], a2[i], 1.e-8);
		BOOST_CHECK_CLOSE(b[i], b2[i], 1.e-8);
	}
	
	Rational_Approximation approx_file(boost::unit_test::framework::master_test_suite().argv[1]);
	int ord3 = approx_file.Get_order();
	hmc_float a03 = approx_file.Get_a0();
	std::vector<hmc_float> a3 = approx_file.Get_a();
	std::vector<hmc_float> b3 = approx_file.Get_b();
	BOOST_CHECK_CLOSE(a0, a03, 1.e-8);
	BOOST_REQUIRE_EQUAL(ord, ord3);
	for(int i=0; i<ord; i++){
		BOOST_CHECK_CLOSE(a[i], a3[i], 1.e-8);
		BOOST_CHECK_CLOSE(b[i], b3[i], 1.e-8);
	}
	
	logger.info() << "Test done!";
}


BOOST_AUTO_TEST_CASE(errors)
{
	using namespace physics::algorithms;
	
	//Calculating 3 different precision rational approximations of x^(-1/2) and comparing errors
	//with those declared by Martin Lüscher in Computational Strategies in Lattice QCD (page 42).
	//Actually, we put here the reference values from the Lüscher's code available at the link
	//http://cern.ch/luscher/ and written to approximate the sign(x) function. But since such
	//a function is there obtained as x*R(x^2), with R(x) the rational approximation of the 
	//inverse square root of x, we can use that code to extract the errors of the approximation.
	hmc_float errors[3];
	errors[0] = 5.0287212083799089e-04;
	errors[1] = 1.2615414383737566e-07;
	errors[2] = 7.9394311806793910e-15;
	//5.e-4, 1.3e-7, 8.e-15};
	Rational_Approximation approx1(6,1,2,1e-5,1);
	//logger.info() << "  Error of the approximation: " << std::setprecision(16) << approx1.Get_error();
	logger.info() << "       Reference value (by Lüscher): " << std::setprecision(16) << errors[0];
	BOOST_CHECK_CLOSE(approx1.Get_error(), errors[0], 1.e-8);
	Rational_Approximation approx2(12,1,2,1e-5,1);
	//logger.info() << "  Error of the approximation: " << std::setprecision(16) << approx2.Get_error();
	logger.info() << "        Reference value (by Lüscher): " << std::setprecision(16) << errors[1];
	BOOST_CHECK_CLOSE(approx2.Get_error(), errors[1], 1.e-8);
	Rational_Approximation approx3(24,1,2,1e-5,1);
	//logger.info() << "  Error of the approximation: " << std::setprecision(16) << approx3.Get_error();
	logger.info() << "        Reference value (by Lüscher): " << std::setprecision(16) << errors[2];
	BOOST_CHECK_CLOSE(approx3.Get_error(), errors[2], 1.e-8);
}


BOOST_AUTO_TEST_CASE(coefficients)
{
	using namespace physics::algorithms;
	
	//Again, we used the Martin Lüscher's code described above to build a test on
	//the coefficients of the rational approximation. Indeed, with such a code we
	//can generate the Approximation of the inverse square root of x in the form
	//(2.97) of the paper Computational Strategies in Lattice QCD, i.e. as ratio
	//of polynomials. To develop the following test we deduced the partial fractions
	//expansion with the Apart function of Mathematica.
	hmc_float A = 0.2801782460422457;
	hmc_float res[5] = {0.0194661374678695, 0.03561358805557574, 0.0821572765515051, 0.21113622523593300, 0.7946025292155642};
	hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145, 0.12517586269872370, 1.0029328743375700};
	hmc_float delta = 5.3847952988591471e-05;
	std::vector<hmc_float> a, b;
	Rational_Approximation approx(5,1,2,0.001,1);
	a = approx.Get_a();
	b = approx.Get_b();
	BOOST_CHECK_CLOSE(approx.Get_a0(), A, 1.e-6);
	for(int i=0; i<5; i++){
	  BOOST_CHECK_CLOSE(a[i], res[i], 1.e-6);
	  BOOST_CHECK_CLOSE(b[i], pol[i], 1.e-6);
	}
	BOOST_CHECK_CLOSE(approx.Get_error(), delta, 1.e-6);
}


BOOST_AUTO_TEST_CASE(rescale)
{
	using namespace physics::algorithms;
	using namespace physics::lattices;
	
	Rational_Approximation approx(15,1,4,1e-5,1,false);
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.567", "--conservative=false", "--num_dev=1"};
	meta::Inputparameters params(6, _params);
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	
	//Min and max eigenvalues for conservative and not conservative case
	hmc_float minEigenvalue = 0.3485318319429664;
	hmc_float maxEigenvalue = 5.2827906935473500;
	hmc_float minEigenvalueCons = 0.321489;
	hmc_float maxEigenvalueCons = 5.546930228224717;

	//Reference rescaled coefficients
	hmc_float a0_ref = 3.78396627036665123;
	hmc_float a_ref[15] = {-2.722986658932525683e-07, -1.2631475639728917521e-06,
	                       -4.360435755914117958e-06, -1.4160691433260296987e-05,
			       -4.5207926211912900696e-05, -0.00014352921651252466598,
			       -0.00045497933178524781715, -0.0014431546521933168083,
			       -0.0045926953908840342788, -0.014747783330565073998,
			       -0.048456946841957317107, -0.16880046472141346792,
			       -0.68431552061715394952, -4.2198332136416603078,
			       -117.03837887995429412};
	hmc_float b_ref[15] = {9.2907369101588763806e-06, 5.0627210699159516975e-05,
	                       0.00016198783558680096995, 0.00044474583457844121164,
			       0.0011563267279465730842, 0.0029447969579766914219,
			       0.007440491639209977949, 0.018751086281155571189,
			       0.047271072978110562079, 0.11959352757092110708,
			       0.30563294009891067704, 0.80189441104810432748,
			       2.2562393220051499831, 7.7746421837770229857,
			       59.252309299420609534};
	//Reference rescaled coefficients conservative
	hmc_float a0_ref_cons = 3.8304052181004228927;
	hmc_float a_ref_cons[15] = {-2.8942286130576286221e-07, -1.3425838169907005556e-06,
				    -4.6346528686679587346e-06, -1.5051222595004849162e-05,
				    -4.8050941840019680289e-05, -0.00015255541700051565286,
				    -0.00048359186633654911497, -0.0015339111096058541221,
				    -0.0048815187425726713766, -0.015675235262161982264,
				    -0.051504285410779379606, -0.17941591204330820108,
				    -0.72735044574402074602, -4.4852081772741634325,
				    -124.39863554566063897};
	hmc_float b_ref_cons[15] = {9.7552862512145591676e-06, 5.3158639325027219212e-05,
				    0.00017008744523117717558, 0.00046698372446689161283,
				    0.0012141446195419417376, 0.003092040766471044512,
				    0.007812526228236587808, 0.019688665814428581158,
				    0.04963469020409881638, 0.12557336479641045823,
				    0.32091499816392937694, 0.84199021010590602287,
				    2.3690543226274733968, 8.1633847494467222106,
				    62.215004455600926292};
	
	Rational_Coefficients coeff = approx.Rescale_Coefficients(minEigenvalue, maxEigenvalue);
	Rational_Coefficients coeff_cons = approx.Rescale_Coefficients(minEigenvalueCons, maxEigenvalueCons);
	
	int ord = coeff.Get_order();
	std::vector<hmc_float> a = coeff.Get_a();
	std::vector<hmc_float> b = coeff.Get_b();
	
	std::vector<hmc_float> a_cons = coeff_cons.Get_a();
	std::vector<hmc_float> b_cons = coeff_cons.Get_b();
	
	//Test result: note that the precision is not so high since
	//the reference code uses a slightly different method to calculate
	//maximum and minimum eigenvalues (I tuned a bit the ref.code adapting the number
	//of loop iterations in finding the max and min eigenvalues, but not too much)
	BOOST_CHECK_CLOSE(coeff.Get_a0(), a0_ref, 5.e-5);
	BOOST_CHECK_CLOSE(coeff_cons.Get_a0(), a0_ref_cons, 5.e-5);
	for(int i=0; i<ord; i++){
		BOOST_CHECK_CLOSE(a[i], a_ref[i], 5.e-5);
		BOOST_CHECK_CLOSE(b[i], b_ref[i], 5.e-5);
		BOOST_CHECK_CLOSE(a_cons[i], a_ref_cons[i], 2.e-4);
		BOOST_CHECK_CLOSE(b_cons[i], b_ref_cons[i], 2.e-4);
	}
	
	logger.info() << "Test done!";
}

BOOST_AUTO_TEST_CASE(input_output)
{
	using namespace physics::algorithms;
	
	int order = 5;
	bool inv = true;         
	int y = 1;            
	int z = 2;            
	int precision = 113;
	hmc_float low = 0.001;    
	hmc_float high = 1.0;
	
	std::vector<hmc_float> a, b;
	hmc_float a0, error;  
	Rational_Approximation approx(order, y, z, low, high, inv, precision);
	a = approx.Get_a();
	b = approx.Get_b();
	a0 = approx.Get_a0();
	error = approx.Get_error();
	
	std::string outputFileName="temporalFileForRationalApproximation";
	BOOST_REQUIRE_MESSAGE(boost::filesystem::exists( outputFileName ) == false,
			      "The file \"" << outputFileName << "\" exists for some reason! Aborting...");
	approx.Save_rational_approximation("temporalFileForRationalApproximation");
	BOOST_REQUIRE_EQUAL(boost::filesystem::exists( outputFileName ), true);
	Rational_Approximation approx2(outputFileName);
	
	BOOST_REQUIRE_EQUAL(order, approx2.Get_order());
	BOOST_REQUIRE_EQUAL(-(double)y/z, approx2.Get_exponent());
	BOOST_REQUIRE_EQUAL(low, approx2.Get_lower_bound());
	BOOST_REQUIRE_EQUAL(high, approx2.Get_upper_bound());
	BOOST_REQUIRE_CLOSE(a0, approx2.Get_a0(), 1.e-13);
	for(int i=0; i<order; i++){
	  BOOST_REQUIRE_CLOSE(a[i], approx2.Get_a()[i], 1.e-13);
	  BOOST_REQUIRE_CLOSE(b[i], approx2.Get_b()[i], 1.e-13);
	}
	BOOST_REQUIRE_CLOSE(error, approx2.Get_error(), 1.e-13);
	
	boost::filesystem::remove(outputFileName);
}
