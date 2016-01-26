/** @file
 * Unit test for the physics::lattices::Spinorfield class
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

#include "spinorfield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Spinorfield
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <cmath>
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
}

BOOST_AUTO_TEST_CASE(fields)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);

	auto fields = create_spinorfields(system, 2, interfacesHandler);

	BOOST_CHECK_EQUAL(fields.size(), 2u);

	release_spinorfields(fields);
}

BOOST_AUTO_TEST_CASE(gamma5)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.zero();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
	sf.cold();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), 1., .1);
	physics::lattices::sax(&sf, { -.5, .3}, sf);
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), .34, .1);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.gaussian(prng);
	sf.zero();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
}

BOOST_AUTO_TEST_CASE(cold)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.gaussian(prng);
	sf.cold();
	BOOST_CHECK_CLOSE(squarenorm(sf), 1., .1);
}


BOOST_AUTO_TEST_CASE(gaussian)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.zero();
	sf.gamma5();
	hmc_float const gamma5 = squarenorm(sf);
	sf.gaussian(prng);
	BOOST_CHECK_NE(squarenorm(sf), gamma5);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	sf.zero();
	sf.gamma5();
	hmc_float const gamma5 = physics::lattices::squarenorm(sf);
	BOOST_REQUIRE_CLOSE(gamma5, 0., .1);
	sf.gaussian(prng);
	BOOST_CHECK_NE(physics::lattices::squarenorm(sf), gamma5);
	sf.zero();
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0., .1);
}

BOOST_AUTO_TEST_CASE(scalar_product)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.gaussian(prng);

	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.zero();

	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();

	Spinorfield gamma(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gamma.zero();
	gamma.gamma5();
	gamma.gamma5();

	const hmc_complex gaussian_scalar_prod = physics::lattices::scalar_product(gaussian, gaussian);
	const hmc_float gaussian_squarenorm = physics::lattices::squarenorm(gaussian);
	BOOST_CHECK_CLOSE(gaussian_scalar_prod.re, gaussian_squarenorm, .1);
	BOOST_CHECK_CLOSE(gaussian_scalar_prod.im, 0., .1);
	const hmc_complex gaussian_scalar_cold = physics::lattices::scalar_product(gaussian, cold);
	const hmc_complex cold_scalar_gaussian = physics::lattices::scalar_product(cold, gaussian);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.re), std::abs(cold_scalar_gaussian.re), .1);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.im), std::abs(cold_scalar_gaussian.im), .1);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(gamma, zero), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(zero, gamma), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(gamma, cold), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(cold, gamma), hmc_complex_zero);
}

BOOST_AUTO_TEST_CASE(sax)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield orig_sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	orig_sf.gaussian(prng);
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::sax(&sf, {.5, 0}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .25 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {2., 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {0., 0.}, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0.);

	orig_sf.cold();
	physics::lattices::sax(&sf, { -.8, .7}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 1.1299999999999968, .1);
	physics::lattices::sax(&sf, {.65, .3}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.51250000000000162, .1);
}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.gaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.zero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::saxpy(&sf, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxpy(&sf, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxpy(&sf, {.3, .1}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .1, .1);
}

BOOST_AUTO_TEST_CASE(saxsbypz)
{
	using physics::lattices::Spinorfield;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Spinorfield gaussian(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	gaussian.gaussian(prng);
	Spinorfield cold(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	cold.cold();
	Spinorfield zero(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	zero.zero();
	Spinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

	physics::lattices::saxsbypz(&sf, {1., 0.}, gaussian, {0., 0.}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, cold, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {.3, .7}, cold, {1., 0.}, zero, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 2.18, .1);
	physics::lattices::saxsbypz(&sf, {.1, .3}, zero, {.7, .3}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 2.98, .1);
}
