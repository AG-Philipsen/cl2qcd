/** @file
 * Implementation of the algorithm to find min and max eigenvalues of an operator
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "find_minmax_eigenvalue.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/util.hpp"
#include "../../logger.hpp"


// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::find_minmax_eigenvalue
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(max)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float ref_max_eig = 5.9887256245527069609;
	
	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, 1.01335);
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	
	hmc_float max = find_max_eigenvalue(matrix, gf, system, 1.e-5);
	
	logger.info() << "ref_max_eig = " << std::setprecision(16) << ref_max_eig;
	logger.info() << "    max_eig = " << std::setprecision(16) << max;
	
	//The precision of this test is not 1.e-8 because the method is not exactly
	//in the same way implemented in the Ref.Code
	BOOST_CHECK_CLOSE(ref_max_eig, max, 1.e-6);
}

BOOST_AUTO_TEST_CASE(min)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float ref_min_eig = 1.053927941164244;
	
	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, 1.01335);
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	
	hmc_float min = find_min_eigenvalue(matrix, gf, system, 1.e-3);
	
	logger.info() << "mass squared = " << std::setprecision(16) << matrix.get_mass()*matrix.get_mass();
	logger.info() << " ref_min_eig = " << std::setprecision(16) << ref_min_eig;
	logger.info() << "     min_eig = " << std::setprecision(16) << min;
	
	BOOST_REQUIRE_SMALL(matrix.get_mass()*matrix.get_mass(), min);
	//The precision of this test is not so high because the method to calculate
	// the eigenvalue is not exactly the same as that implemented in the Ref.Code
	BOOST_CHECK_CLOSE(ref_min_eig, min, 1.e-3);
}