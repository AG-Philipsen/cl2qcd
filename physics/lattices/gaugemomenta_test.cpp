/** @file
 * Unit test for the physics::lattices::Gaugemomenta class
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "gaugemomenta.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugemomenta
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include "../../hardware/device.hpp"
#include "../../hardware/code/gaugemomentum.hpp"
#include "util.hpp"
#include "../../meta/util.hpp"
#include "../../interfaceImplementations/latticesParameters.hpp"
#include "../../interfaceImplementations/physicsParameters.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

static void fill_buffer(const hardware::buffers::Gaugemomentum * buf, int seed);

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system, gaugemomentaParametersInterface);

	BOOST_REQUIRE_NE(gm.get_buffers().size(), 0u);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system, gaugemomentaParametersInterface);

	// fill all buffers with noise
	for(auto buffer: gm.get_buffers()) {
		fill_buffer(buffer, 13);
	}

	// run code
	gm.zero();

	// check result
	for(auto buffer: gm.get_buffers()) {
		size_t num_elems = buffer->get_elements();
		ae * host_mem = new ae[num_elems];
		buffer->get_device()->getGaugemomentumCode()->exportGaugemomentumBuffer(host_mem, buffer);
		const ae zero = {0, 0, 0, 0, 0, 0, 0, 0};
		for(size_t i = 0; i < num_elems; ++i) {
			BOOST_REQUIRE_EQUAL(host_mem[i], zero);
		}
	}
}

static void fill_buffer(const hardware::buffers::Gaugemomentum * buf, int seed)
{
	size_t num_elems = buf->get_elements();
	ae * host_mem = new ae[num_elems];
	fill(host_mem, num_elems, seed);
	buf->get_device()->getGaugemomentumCode()->importGaugemomentumBuffer(buf, host_mem);
}

BOOST_AUTO_TEST_CASE(gaussian)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system, gaugemomentaParametersInterface);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	// fill with zeros
	gm.zero();

	// run code
	gm.gaussian(prng);

	// simple verification
	BOOST_REQUIRE_NE(squarenorm(gm), 0.);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using namespace physics::lattices;
	using physics::lattices::squarenorm;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system, gaugemomentaParametersInterface);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	// only two very simple tests

	// this must be zero...
	gm.zero();
	BOOST_CHECK_EQUAL(squarenorm(gm), 0.);

	// this should never be zero
	gm.gaussian(prng);
	BOOST_CHECK_NE(squarenorm(gm), 0.);

	// and the same for the asynchroneous variant
	Scalar<hmc_float> res(system);

	gm.zero();
	squarenorm(&res, gm);
	BOOST_CHECK_EQUAL(res.get(), 0.);

	// this should never be zero
	gm.gaussian(prng);
	squarenorm(&res, gm);
	BOOST_CHECK_NE(res.get(), 0.);

	pseudo_randomize<Gaugemomenta, ae>(&gm, 5);
	squarenorm(&res, gm);
	BOOST_CHECK_CLOSE(res.get(), 5509.4389078650529, .01);

	pseudo_randomize<Gaugemomenta, ae>(&gm, 51);
	squarenorm(&res, gm);
	BOOST_CHECK_CLOSE(res.get(), 5484.798507726874, .01);
}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Gaugemomenta;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Gaugemomenta gauss(system, gaugemomentaParametersInterface);
	gauss.gaussian(prng);
	Gaugemomenta zz(system, gaugemomentaParametersInterface);
	zz.zero();
	Gaugemomenta gm(system, gaugemomentaParametersInterface);

	physics::lattices::saxpy(&gm, 1., gauss, zz);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(gm), physics::lattices::squarenorm(gauss), 1.e-8);
	physics::lattices::saxpy(&gm, 0., gauss, zz);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(gm), 0);
}

BOOST_AUTO_TEST_CASE(halo_update)
{
	using namespace physics::lattices;

	hmc_float orig_squarenorm, new_squarenorm;

	// simple test, squarenorm should not get changed by halo exchange
	const char * _params[] = {"foo", "--ntime=16"};
	meta::Inputparameters params(2, _params);
	physics::lattices::GaugemomentaParametersImplementation gaugemomentaParametersInterface{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	const Gaugemomenta gm(system, gaugemomentaParametersInterface);

	gm.gaussian(prng);
	orig_squarenorm = physics::lattices::squarenorm(gm);
	gm.update_halo();
	new_squarenorm = physics::lattices::squarenorm(gm);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	gm.zero();
	orig_squarenorm = physics::lattices::squarenorm(gm);
	gm.update_halo();
	new_squarenorm = physics::lattices::squarenorm(gm);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	gm.gaussian(prng);
	orig_squarenorm = physics::lattices::squarenorm(gm);
	gm.update_halo();
	new_squarenorm = physics::lattices::squarenorm(gm);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	pseudo_randomize<Gaugemomenta, ae>(&gm, 5);
	orig_squarenorm = physics::lattices::squarenorm(gm);
	gm.update_halo();
	new_squarenorm = physics::lattices::squarenorm(gm);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	pseudo_randomize<Gaugemomenta, ae>(&gm, 51);
	orig_squarenorm = physics::lattices::squarenorm(gm);
	gm.update_halo();
	new_squarenorm = physics::lattices::squarenorm(gm);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);
}
