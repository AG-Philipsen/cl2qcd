/** @file
 * Unit test for the physics::lattices::Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "staggeredfield_eo.hpp"
#include "util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Staggeredfield_eo
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <cmath>


BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Staggeredfield_eo sf(system);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo sf(system);
	sf.set_zero();
	hmc_float const sq = physics::lattices::squarenorm(sf);
	BOOST_REQUIRE_EQUAL(sq, 0);
	sf.set_gaussian(prng);
	BOOST_CHECK_NE(physics::lattices::squarenorm(sf), sq);
	sf.set_zero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo sf(system);
	sf.set_gaussian(prng);
	sf.set_zero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
}

BOOST_AUTO_TEST_CASE(cold)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo sf(system);
	sf.set_gaussian(prng);
	sf.set_cold();
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.5, 1.e-8);
}

BOOST_AUTO_TEST_CASE(gaussian)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo sf(system);
	sf.set_cold();
	hmc_float const sq = physics::lattices::squarenorm(sf);
	sf.set_gaussian(prng);
	BOOST_CHECK_NE(physics::lattices::squarenorm(sf), sq);
}

BOOST_AUTO_TEST_CASE(scalar_product)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo gaussian(system);
	gaussian.set_gaussian(prng);

	Staggeredfield_eo zero(system);
	zero.set_zero();

	Staggeredfield_eo cold(system);
	cold.set_cold();

	const hmc_complex gaussian_scalar_prod = physics::lattices::scalar_product(gaussian, gaussian);
	const hmc_float gaussian_squarenorm = physics::lattices::squarenorm(gaussian);
	BOOST_CHECK_CLOSE(gaussian_scalar_prod.re, gaussian_squarenorm, 1.e-8);
	BOOST_CHECK_EQUAL(gaussian_scalar_prod.im, 0);
	const hmc_complex gaussian_scalar_cold = physics::lattices::scalar_product(gaussian, cold);
	const hmc_complex cold_scalar_gaussian = physics::lattices::scalar_product(cold, gaussian);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.re), std::abs(cold_scalar_gaussian.re), 1.e-8);
	BOOST_CHECK_CLOSE(std::abs(gaussian_scalar_cold.im), std::abs(cold_scalar_gaussian.im), 1.e-8);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(gaussian, zero), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(zero, gaussian), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(zero, cold), hmc_complex_zero);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product(cold, zero), hmc_complex_zero);
}

BOOST_AUTO_TEST_CASE(sax)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo orig_sf(system);
	orig_sf.set_gaussian(prng);
	Staggeredfield_eo sf(system);

	physics::lattices::sax(&sf, {0.5, 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25 * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, {2., 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4. * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, {0., 0.}, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);

	orig_sf.set_cold();
	physics::lattices::sax(&sf, { -.8, .7}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.565, 1.e-8);
	physics::lattices::sax(&sf, {.65, .3}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25625, 1.e-8);

}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo gaussian(system);
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system);
	cold.set_cold();
	Staggeredfield_eo zero(system);
	zero.set_zero();
	Staggeredfield_eo sf(system);

	physics::lattices::saxpy(&sf, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpy(&sf, {0., 0.}, gaussian, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(cold), 1.e-8);
	physics::lattices::saxpy(&sf, {.3, .1}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .85, 1.e-8);
}

BOOST_AUTO_TEST_CASE(saxpby)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo gaussian(system);
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system);
	cold.set_cold();
	Staggeredfield_eo zero(system);
	zero.set_zero();
	Staggeredfield_eo sf(system);

	physics::lattices::saxpby(&sf, {1., 0.}, gaussian, {0., 0.}, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, cold, {1., 0.}, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	physics::lattices::saxpby(&sf, {.3, .7}, cold, {1., 0.}, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.29, 1.e-8);
	physics::lattices::saxpby(&sf, {.1, .3}, zero, {.56, .65}, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.36805, 1.e-8);
}

BOOST_AUTO_TEST_CASE(saxpbypz)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Staggeredfield_eo gaussian(system);
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system);
	cold.set_cold();
	Staggeredfield_eo zero(system);
	zero.set_zero();
	Staggeredfield_eo sf(system);

	physics::lattices::saxpbypz(&sf, {1., 0.}, gaussian, {0., 0.}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpbypz(&sf, {0., 0.}, cold, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpbypz(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpbypz(&sf, {.3, .7}, cold, {1., 0.}, zero, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 1.09, 1.e-8);
	physics::lattices::saxpbypz(&sf, {.1, .3}, zero, {.7, .3}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 1.49, 1.e-8);
}

BOOST_AUTO_TEST_CASE(pseudorandomize)
{
	using physics::lattices::Staggeredfield_eo;
	
	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	Staggeredfield_eo sf(system);
	sf.set_zero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	physics::lattices::pseudo_randomize<Staggeredfield_eo, su3vec>(&sf, 123);
	logger.info() << "The squarenorm of the pseudorandomized field is " << physics::lattices::squarenorm(sf);
}

//To be added (so far only single GPU and only EO preconditioning)...
#if 0

BOOST_AUTO_TEST_CASE(conversion)
{
	using physics::lattices::Spinorfield;
	using physics::lattices::Spinorfield_eo;
	using physics::lattices::squarenorm;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	for(size_t i = 0; i < 11; i++) {
		const Spinorfield orig(system);
		orig.gaussian(prng);
		const Spinorfield recreated(system);
		recreated.gaussian(prng);
		const Spinorfield_eo even(system);
		const Spinorfield_eo odd(system);
		even.gaussian(prng);
		odd.gaussian(prng);

		convert_to_eoprec(&even, &odd, orig);
		log_squarenorm("even: ", even);
		log_squarenorm("odd: ", odd);
		log_squarenorm("orig: ", orig);
		BOOST_CHECK_CLOSE(squarenorm(even) + squarenorm(odd), squarenorm(orig), .1);

		convert_from_eoprec(&recreated, even, odd);
		log_squarenorm("recreated: ", recreated);
		BOOST_CHECK_CLOSE(squarenorm(recreated), squarenorm(orig), .1);
	}
}

BOOST_AUTO_TEST_CASE(halo_update)
{
	using namespace physics::lattices;

	hmc_float orig_squarenorm, new_squarenorm;

	// simple test, squarenorm should not get changed by halo exchange
	const char * _params[] = {"foo", "--ntime=16"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	const Spinorfield_eo sf(system);

	sf.gaussian(prng);
	orig_squarenorm = physics::lattices::squarenorm(sf);
	sf.update_halo();
	new_squarenorm = physics::lattices::squarenorm(sf);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	sf.zero();
	orig_squarenorm = physics::lattices::squarenorm(sf);
	sf.update_halo();
	new_squarenorm = physics::lattices::squarenorm(sf);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);

	sf.cold();
	orig_squarenorm = physics::lattices::squarenorm(sf);
	sf.update_halo();
	new_squarenorm = physics::lattices::squarenorm(sf);
	BOOST_CHECK_EQUAL(orig_squarenorm, new_squarenorm);
}

#endif
