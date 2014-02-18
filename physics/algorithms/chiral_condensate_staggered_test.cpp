/** @file
 * Tests of the flavour doublets algorithms
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "chiral_condensate_staggered.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::chiral_condensate_stagg
#include <boost/test/unit_test.hpp>

#include "../prng.hpp"
#include "../../host_functionality/logger.hpp"

void test_chiral_condensate_stagg(std::string content, hmc_complex pbp_ref, bool cold, bool thetaT)
{
	using namespace physics::lattices;

	std::vector<const char*> options(1, "foo");
	options.push_back("--ntime=4");
	options.push_back("--fermact=rooted_stagg");
	options.push_back("--num_tastes=2");
	options.push_back("--mass=0.1");
	options.push_back("--theta_fermion_spatial=0");
	if(thetaT)
	  options.push_back("--theta_fermion_temporal=1");
	else
	  options.push_back("--theta_fermion_temporal=0");
	options.push_back("--num_sources=10");
 	options.push_back("--sourcetype=volume");
	std::string tmp = "--sourcecontent=" + content;
	options.push_back(tmp.c_str());
	
	meta::Inputparameters params(10, &(options[0]));
	hardware::System system(params);
	physics::PRNG prng(system);
	const Gaugefield *gf;
	if(cold)
	  gf = new Gaugefield(system, prng, false);
	else //This configuration for the Ref.Code is the same as for example dks_input_5
	  gf = new Gaugefield(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	
	hmc_complex pbp = physics::algorithms::chiral_condensate_staggered(*gf, prng, system);
	
	logger.info() << "Chiral condensate pbp = (" << std::setprecision(12) << pbp.re << "," << pbp.im << ")";
	BOOST_CHECK_CLOSE(pbp.re, pbp_ref.re, 1.e-8);
	BOOST_CHECK_CLOSE(pbp.im, pbp_ref.im, 1.e-8);
}


BOOST_AUTO_TEST_CASE(volume_source_one)
{
	test_chiral_condensate_stagg("one", {15.0, 0.0}, true, 0);
	test_chiral_condensate_stagg("one", {0.29411764705882464943, 0.0}, true, 1);
// 	test_chiral_condensate_stagg("one", {1.0,0.0}, false, 0);
// 	test_chiral_condensate_stagg("one", {1.0,0.0}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_z4)
{
// 	test_chiral_condensate_stagg("z4", {1.0,0.0}, true, 0);
// 	test_chiral_condensate_stagg("z4", {1.0,0.0}, true, 1);
// 	test_chiral_condensate_stagg("z4", {1.0,0.0}, false, 0);
// 	test_chiral_condensate_stagg("z4", {1.0,0.0}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_gaussian)
{
// 	test_chiral_condensate_stagg("gaussian", {1.0,0.0}, true, 0);
// 	test_chiral_condensate_stagg("gaussian", {1.0,0.0}, true, 1);
// 	test_chiral_condensate_stagg("gaussian", {1.0,0.0}, false, 0);
// 	test_chiral_condensate_stagg("gaussian", {1.0,0.0}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_z2)
{
// 	test_chiral_condensate_stagg("z2", {1.0,0.0}, true, 0);
// 	test_chiral_condensate_stagg("z2", {1.0,0.0}, true, 1);
// 	test_chiral_condensate_stagg("z2", {1.0,0.0}, false, 0);
// 	test_chiral_condensate_stagg("z2", {1.0,0.0}, false, 1);
}
