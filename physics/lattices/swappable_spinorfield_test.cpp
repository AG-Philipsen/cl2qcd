/** @file
 * Tests of the physics::lattices::SwappableSpinorfield class
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

#include "swappable_spinorfield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::SwappableSpinorfield
#include <boost/test/unit_test.hpp>

#include "../../interfaceImplementations/interfacesHandler.hpp"
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
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	logger.debug() << "Devices: " << system.get_devices().size();

	SwappableSpinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
}

BOOST_AUTO_TEST_CASE(swap)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	logger.debug() << "Devices: " << system.get_devices().size();

	SwappableSpinorfield sf(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
	BOOST_REQUIRE_LT(0, sf.get_buffers().size());

	sf.swap_out();
	BOOST_REQUIRE_EQUAL(0, sf.get_buffers().size());

	sf.swap_in();
	BOOST_REQUIRE_LT(0, sf.get_buffers().size());
}

