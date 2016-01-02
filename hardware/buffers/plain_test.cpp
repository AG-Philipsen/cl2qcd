/** @file
 * Testcases for the hardware::buffers::Plain template
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

#include "plain.hpp"
#include "../system.hpp"
#include "../device.hpp"
#include "../../common_header_files/types.h"
#include "../../meta/type_ops.hpp"
#include "../interfaceMockups.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Plain
#include <boost/test/unit_test.hpp>

template<typename T> void test(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;

	Plain<T> dummy(elems, device);
	BOOST_REQUIRE_EQUAL(dummy.get_elements(), elems);
	BOOST_REQUIRE_EQUAL(dummy.get_bytes(), elems * sizeof(T));
	const cl_mem * tmp = dummy;
	BOOST_REQUIRE(tmp);
	BOOST_REQUIRE(*tmp);
	BOOST_REQUIRE_EQUAL(device->get_id(), dummy.get_device()->get_id());

	T* in = new T[elems];
	T* out = new T[elems];
	fill(in, elems, 1);
	fill(out, elems, 2);

	dummy.load(in);
	dummy.dump(out);
	BOOST_CHECK_EQUAL_COLLECTIONS(in, in + elems, out, out + elems);

	Plain<T> dummy2(elems, device);
	fill(in, elems, 3);
	fill(out, elems, 4);
	dummy.load(in);
	copyData(&dummy2, &dummy);
	dummy2.dump(out);
	BOOST_CHECK_EQUAL_COLLECTIONS(in, in + elems, out, out + elems);

	fill(out, elems, 5);
	hardware::SynchronizationEvent event = dummy2.dump_async(out);
	event.wait();
	BOOST_CHECK_EQUAL_COLLECTIONS(in, in + elems, out, out + elems);

	delete[] in;
	delete[] out;
}

template<typename T> void test(bool requireDouble = false)
{
	using namespace hardware;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const std::vector<Device*>& devices = system.get_devices();
	for(Device * device : devices)
	{
		if(!requireDouble || device->is_double_supported()) {
			test<T>(1, device);
			test<T>(1024, device);
		}
	}
}

BOOST_AUTO_TEST_CASE(float_buffer)
{
	test<float>();
}

BOOST_AUTO_TEST_CASE(double_buffer)
{
	test<double>(true);
}

BOOST_AUTO_TEST_CASE(hmc_complex_buffer)
{
	test<hmc_complex>(true);
}

BOOST_AUTO_TEST_CASE(SU3_buffer)
{
	test<Matrixsu3>(true);
}
