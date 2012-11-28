/** @file
 * Tests of the hardware::SynchronizationEvent class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::SynchronizationEvent
#include <boost/test/unit_test.hpp>

#include "synchronization_event.hpp"
#include "system.hpp"

BOOST_AUTO_TEST_CASE(good_case)
{
	using hardware::SynchronizationEvent;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);

	cl_int err;
	cl_event event = clCreateUserEvent(system, &err);
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
}

BOOST_AUTO_TEST_CASE(error_case)
{
	using hardware::SynchronizationEvent;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);

	cl_int err;
	cl_event event = clCreateUserEvent(system, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	SynchronizationEvent se1(event);
	BOOST_REQUIRE_EQUAL(se1.is_finished(), false);

	err = clSetUserEventStatus(event, -1);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	BOOST_CHECK_THROW(se1.is_finished(), hardware::OpenclException);
	BOOST_CHECK_THROW(se1.wait(), hardware::OpenclException);

	SynchronizationEvent se3;
	BOOST_CHECK_THROW(se3.is_finished(), hardware::OpenclException);
	BOOST_CHECK_THROW(se3.wait(), hardware::OpenclException);
}
