/** @file
 * Unit test for the physics::lattices::Gaugemomenta class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomenta.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugemomenta
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"
#include "../../meta/type_ops.hpp"

static void fill_buffer(const hardware::buffers::Gaugemomentum * buf, int seed);

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system);

	BOOST_REQUIRE_NE(gm.get_buffers().size(), 0u);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system);

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
		buffer->get_device()->get_gaugemomentum_code()->exportGaugemomentumBuffer(host_mem, buffer);
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
	buf->get_device()->get_gaugemomentum_code()->importGaugemomentumBuffer(buf, host_mem);
}

BOOST_AUTO_TEST_CASE(gaussian)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gm(system);
	physics::PRNG prng(system);

	// fill with zeros
	gm.zero();

	// run code
	gm.gaussian(prng);

	BOOST_FAIL("Not implemented");
}

