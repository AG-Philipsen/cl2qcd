/** @file
 * Tests of the multi-shifted inverter algorithm
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "solver_shifted.hpp"
#include "../../logger.hpp"
#include "../lattices/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::solvers::solver_shifted
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(cgm)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	
	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//This are some possible values of sigma
	hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145,
	                                        0.12517586269872370, 1.0029328743375700};
	std::vector<hmc_float> sigma(pol, pol + sizeof(pol)/sizeof(hmc_float));
	physics::fermionmatrix::MdagM_eo matrix(system, 0.567);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 13);
	
	int iter=cg_m(out, sigma, matrix, gf, b, system, 1.e-8);
}



