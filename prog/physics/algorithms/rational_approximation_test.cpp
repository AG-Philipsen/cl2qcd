/** @file
 * Tests of the Rational_Approximation algorithm
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "rational_approximation.hpp"
#include "../../logger.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::rational_approximation
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(initialization)
{
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
	hmc_float res[5] = {0.0194661374678695, 0.03561358805557574, 0.0821572765515051,
	                                        0.21113622523593300, 0.7946025292155642};
	hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145,
	                                        0.12517586269872370, 1.0029328743375700};
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
	logger.info() << approx;
	
	getchar();
	
	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, 1.01335);
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	
	Rational_Coefficients coeff = approx.Rescale_Coefficients(matrix, gf, system, 1.e-3);
	
	int ord2 = coeff.Get_order();
	hmc_float a02 = coeff.Get_a0();
	std::vector<hmc_float> a2 = coeff.Get_a();
	std::vector<hmc_float> b2 = coeff.Get_b();
	
	logger.info() << "a02 = " << a02;
	for(int i=0; i<ord2; i++){
		logger.info() << "a2[" << i << "] = " << a2[i];
		logger.info() << "b2[" << i << "] = " << b2[i];
	}
	
	logger.info() << "Test done!";
}


