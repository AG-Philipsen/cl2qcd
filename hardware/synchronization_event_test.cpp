/** @file
 * Tests of the hardware::SynchronizationEvent class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::SynchronizationEvent
#include <boost/test/unit_test.hpp>

#include "synchronization_event.hpp"
#include "system.hpp"
#include "interfaceMockups.hpp"

BOOST_AUTO_TEST_CASE(good_case)
{
	using hardware::SynchronizationEvent;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );

	cl_int err;
	cl_event event = clCreateUserEvent(system.getContext(), &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	SynchronizationEvent se1(event);
	BOOST_CHECK_EQUAL(se1.is_finished(), false);

	SynchronizationEvent se2(se1);
	BOOST_CHECK_EQUAL(se1.is_finished(), false);
	BOOST_CHECK_EQUAL(se2.is_finished(), false);

	err = clSetUserEventStatus(event, CL_COMPLETE);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	BOOST_CHECK_EQUAL(se1.is_finished(), true);
	BOOST_CHECK_EQUAL(se2.is_finished(), true);

	err = clReleaseEvent(event);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	BOOST_CHECK_EQUAL(se1.is_finished(), true);
	BOOST_CHECK_EQUAL(se2.is_finished(), true);

	se1 = se2;

	BOOST_CHECK_EQUAL(se1.is_finished(), true);
	BOOST_CHECK_EQUAL(se2.is_finished(), true);

	// this should not stall
	se1.wait();
	se2.wait();

	SynchronizationEvent se3;
	SynchronizationEvent se4 = se3;
	se3 = se4;

	auto const raw_events = hardware::get_raw_events({se1, se3});
	BOOST_CHECK_EQUAL(raw_events.size(), 1);
	BOOST_CHECK_EQUAL(raw_events[0], event);
}

BOOST_AUTO_TEST_CASE(error_case)
{
	using hardware::SynchronizationEvent;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );

	cl_int err;
	cl_event event = clCreateUserEvent(system.getContext(), &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	SynchronizationEvent se1(event);
	BOOST_REQUIRE_EQUAL(se1.is_finished(), false);

	err = clSetUserEventStatus(event, -1);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	BOOST_CHECK_THROW(se1.is_finished(), hardware::OpenclException);
	BOOST_CHECK_THROW(se1.wait(), hardware::OpenclException);

	SynchronizationEvent se3;
	BOOST_CHECK_EQUAL(se3.is_finished(), true);
	se3.wait();
}
