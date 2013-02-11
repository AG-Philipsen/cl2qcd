/** @file
 * Unit test for the physics::lattices::Spinorfield_eo class
 */

#include "spinorfield_eo.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Spinorfield_eo
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"
#include "../../meta/type_ops.hpp"
#include <cmath>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Spinorfield_eo sf(system);
}

BOOST_AUTO_TEST_CASE(gamma5)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
	sf.zero();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
	sf.cold();
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), .5, .1);
	physics::lattices::sax(&sf, { -.5, .3}, sf);
	sf.gamma5();
	BOOST_CHECK_CLOSE(squarenorm(sf), .17, .1);
}

BOOST_AUTO_TEST_CASE(zero)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
	sf.gaussian(prng);
	sf.zero();
	BOOST_CHECK_CLOSE(squarenorm(sf), 0., .1);
}

BOOST_AUTO_TEST_CASE(cold)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
	sf.gaussian(prng);
	sf.cold();
	BOOST_CHECK_CLOSE(squarenorm(sf), .5, .1);
}


BOOST_AUTO_TEST_CASE(gaussian)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
	sf.zero();
	sf.gamma5();
	hmc_float const gamma5 = squarenorm(sf);
	sf.gaussian(prng);
	BOOST_CHECK_NE(squarenorm(sf), gamma5);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
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
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo gaussian(system);
	gaussian.gaussian(prng);

	Spinorfield_eo zero(system);
	zero.zero();

	Spinorfield_eo cold(system);
	cold.cold();

	Spinorfield_eo gamma(system);
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
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo orig_sf(system);
	orig_sf.gaussian(prng);
	Spinorfield_eo sf(system);

	physics::lattices::sax(&sf, {.5, 0}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .25 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {2., 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4 * physics::lattices::squarenorm(orig_sf), .1);
	physics::lattices::sax(&sf, {0., 0.}, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0.);

	orig_sf.cold();
	physics::lattices::sax(&sf, { -.8, .7}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.56499999999999906, .1);
	physics::lattices::sax(&sf, {.65, .3}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25625000000000059, .1);

}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo gaussian(system);
	gaussian.gaussian(prng);
	Spinorfield_eo cold(system);
	cold.gaussian(prng);
	Spinorfield_eo zero(system);
	zero.zero();
	Spinorfield_eo sf(system);

	physics::lattices::saxpy(&sf, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxpy(&sf, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxpy(&sf, {.3, .1}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 1543.1431463138993, .1);
}

BOOST_AUTO_TEST_CASE(saxsbypz)
{
	using physics::lattices::Spinorfield_eo;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo gaussian(system);
	gaussian.gaussian(prng);
	Spinorfield_eo cold(system);
	cold.gaussian(prng);
	Spinorfield_eo zero(system);
	zero.zero();
	Spinorfield_eo sf(system);

	physics::lattices::saxsbypz(&sf, {1., 0.}, gaussian, {0., 0.}, cold, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, cold, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), .1);
	physics::lattices::saxsbypz(&sf, {.3, .7}, cold, {1., 0.}, zero, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 6728.1041179285985, .1);
	physics::lattices::saxsbypz(&sf, {.1, .3}, zero, {.7, .3}, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 9197.1331520308358, .1);
}

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
		BOOST_CHECK_CLOSE(squarenorm(even) + squarenorm(odd), squarenorm(orig), .1);

		convert_from_eoprec(&recreated, even, odd);
		BOOST_CHECK_CLOSE(squarenorm(recreated), squarenorm(orig), .1);
	}
}
