/** @file
 * Unit test for the physics::lattices::Scalar<hmc_complex> class
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

#include "vector.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"


// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Scalar<hmc_complex>
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(vector_float)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Vector<hmc_float> foo(3, system);
	std::vector<hmc_float> ref(3, 13.);
	foo.store(ref);
	std::vector<hmc_float> got = foo.get();
	
	BOOST_REQUIRE_EQUAL(foo.get_vector_size(), 3);
	BOOST_REQUIRE_EQUAL(ref.size(), got.size());
	for(uint i=0; i<got.size(); i++)
	  BOOST_REQUIRE_EQUAL(got[i], ref[i]);
}


BOOST_AUTO_TEST_CASE(vector_complex)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Vector<hmc_complex> foo(3, system);
	std::vector<hmc_complex> ref(3, {13., 3.1415});
	foo.store(ref);
	std::vector<hmc_complex> got = foo.get();
	
	BOOST_REQUIRE_EQUAL(foo.get_vector_size(), 3);
	BOOST_REQUIRE_EQUAL(ref.size(), got.size());
	for(uint i=0; i<got.size(); i++){
	  BOOST_REQUIRE_EQUAL(got[i].re, ref[i].re);
	  BOOST_REQUIRE_EQUAL(got[i].im, ref[i].im);
	}
}


BOOST_AUTO_TEST_CASE(vector_char)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Vector<char> foo(3, system);
	std::vector<char> ref(3, 'a');
	foo.store(ref);
	std::vector<char> got = foo.get();
	
	BOOST_REQUIRE_EQUAL(foo.get_vector_size(), 3);
	BOOST_REQUIRE_EQUAL(ref.size(), got.size());
	for(uint i=0; i<got.size(); i++){
	  BOOST_REQUIRE_EQUAL(got[i], ref[i]);
	}
}
