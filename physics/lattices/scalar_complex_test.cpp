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

#include "scalar_complex.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Scalar<hmc_complex>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> foo(system);
	hmc_complex ref = {42., 13.};
	foo.store(ref);

	BOOST_REQUIRE_EQUAL(ref, foo.get());
}

BOOST_AUTO_TEST_CASE(multiplication)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> left(system);
	Scalar<hmc_complex> right(system);
	Scalar<hmc_complex> res(system);

	left.store( {1., 2.});
	right.store( {3., 4.});

	multiply(&res, left, right);

	hmc_complex ref = { -5., 10.};
	BOOST_REQUIRE_EQUAL(ref, res.get());
}

BOOST_AUTO_TEST_CASE(division)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> left(system);
	Scalar<hmc_complex> right(system);
	Scalar<hmc_complex> res(system);

	left.store( {1., 2.});
	right.store( {3., 4.});

	divide(&res, left, right);

	hmc_complex res_host = res.get();
	BOOST_REQUIRE_CLOSE(0.44, res_host.re, 1);
	BOOST_REQUIRE_CLOSE(0.08, res_host.im, 1);
}

BOOST_AUTO_TEST_CASE(conversion)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_float> real_dev(system);
	Scalar<hmc_complex> complex_dev(system);

	hmc_float real_host = 1.37;
	hmc_complex complex_host = {2.34, 5.67};

	real_dev.store(real_host);
	complex_dev.store(complex_host);
	convert(&complex_dev, real_dev);
	hmc_complex ref = {real_host, 0.};
	BOOST_REQUIRE_EQUAL(ref, complex_dev.get());
}
