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

#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <cmath>


BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--fermact=rooted_stagg"};
	meta::Inputparameters params(2, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	logger.debug() << "Devices: " << system.get_devices().size();

	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
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

    auto fields = create_staggeredfields_eo(system, 2, interfacesHandler);

    BOOST_CHECK_EQUAL(fields.size(), 2u);

    release_staggeredfields_eo(fields);
}

BOOST_AUTO_TEST_CASE(squarenorm)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
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

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	sf.set_gaussian(prng);
	sf.set_zero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
}

BOOST_AUTO_TEST_CASE(cold)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	sf.set_gaussian(prng);
	sf.set_cold();
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.5, 1.e-8);
}

BOOST_AUTO_TEST_CASE(gaussian)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	sf.set_cold();
	hmc_float const sq = physics::lattices::squarenorm(sf);
	sf.set_gaussian(prng);
	BOOST_CHECK_NE(physics::lattices::squarenorm(sf), sq);
}

BOOST_AUTO_TEST_CASE(scalar_product)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo gaussian(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	gaussian.set_gaussian(prng);

	Staggeredfield_eo zero(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	zero.set_zero();

	Staggeredfield_eo cold(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
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
	//Real part function
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(gaussian, gaussian),
			   physics::lattices::squarenorm(gaussian));
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(gaussian, cold), gaussian_scalar_cold.re);
	BOOST_CHECK_EQUAL(physics::lattices::scalar_product_real_part(cold, gaussian), cold_scalar_gaussian.re);
}

BOOST_AUTO_TEST_CASE(sax)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo orig_sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	physics::lattices::Scalar<hmc_complex> cplx(system);
	physics::lattices::Scalar<hmc_float> real(system);
	physics::lattices::Vector<hmc_float> real_vec(5, system);
	
	//Complex
	orig_sf.set_gaussian(prng);
	physics::lattices::sax(&sf, {0.5, 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25 * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, {2., 0.}, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4. * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, {0., 0.}, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	orig_sf.set_cold();
	cplx.store({-0.8, 0.7});
	physics::lattices::sax(&sf, cplx, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.565, 1.e-8);
	cplx.store({0.65, 0.3});
	physics::lattices::sax(&sf, cplx, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25625, 1.e-8);
	//Real
	orig_sf.set_gaussian(prng);
	physics::lattices::sax(&sf, 0.5, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.25 * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, 2.0, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 4. * physics::lattices::squarenorm(orig_sf), 1.e-8);
	physics::lattices::sax(&sf, 0.0, orig_sf);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	orig_sf.set_cold();
	real.store(-0.123);
	physics::lattices::sax(&sf, real, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.0075645, 1.e-8);
	real.store(0.653);
	physics::lattices::sax(&sf, real, orig_sf);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.2132045, 1.e-8);
	//Vector
	real_vec.store(std::vector<hmc_float>(5, 0.31415));
	for(int i=0; i<5; i++){
	  physics::lattices::sax(&sf, real_vec, i, orig_sf);
	  BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.04934511125, 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(saxpy)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo gaussian(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	cold.set_cold();
	Staggeredfield_eo zero(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	zero.set_zero();
	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	physics::lattices::Scalar<hmc_complex> cplx(system);
	physics::lattices::Scalar<hmc_float> real(system);
	physics::lattices::Vector<hmc_float> real_vec(5, system);
	cplx.store({0.3, 0.1});
	real.store(0.3);

	//Complex
	physics::lattices::saxpy(&sf, {1., 0.}, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpy(&sf, {0., 0.}, gaussian, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(cold), 1.e-8);
	physics::lattices::saxpy(&sf, cplx, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .85, 1.e-8);
	//Real
	physics::lattices::saxpy(&sf, 1.0, gaussian, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpy(&sf, 0.0, gaussian, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(cold), 1.e-8);
	physics::lattices::saxpy(&sf, real, cold, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), .845, 1.e-8);
	//Vector
	real_vec.store(std::vector<hmc_float>(5, 0.31415));
	for(int i=0; i<5; i++){
	  physics::lattices::saxpy(&sf, real_vec, i, cold, zero);
	  BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.04934511125, 1.e-8);
	  physics::lattices::saxpy(&sf, real_vec, i, cold, cold);
	  BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.86349511125, 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(saxpby)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo gaussian(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	cold.set_cold();
	Staggeredfield_eo zero(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	zero.set_zero();
	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	physics::lattices::Scalar<hmc_complex> cplx(system);
	physics::lattices::Scalar<hmc_float> real(system);
	physics::lattices::Vector<hmc_float> real_vec1(5, system);
	physics::lattices::Vector<hmc_float> real_vec2(5, system);

	//Complex
	physics::lattices::saxpby(&sf, {1., 0.}, gaussian, {0., 0.}, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, cold, {1., 0.}, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	cplx.store({0.3, 0.7});
	physics::lattices::saxpby(&sf, cplx, cold, cplx, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.29, 1.e-8);
	cplx.store({0.56, 0.65});
	physics::lattices::saxpby(&sf, cplx, zero, cplx, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.36805, 1.e-8);
	//Real
	physics::lattices::saxpby(&sf, {1., 0.}, gaussian, {0., 0.}, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, cold, {1., 0.}, gaussian);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), physics::lattices::squarenorm(gaussian), 1.e-8);
	physics::lattices::saxpby(&sf, {0., 0.}, gaussian, {0., 0.}, gaussian);
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	real.store(0.3);
	physics::lattices::saxpby(&sf, real, cold, real, zero);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.045, 1.e-8);
	real.store(-0.56);
	physics::lattices::saxpby(&sf, real, zero, real, cold);
	BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.1568, 1.e-8);
	//Vector
	real_vec1.store(std::vector<hmc_float>(5, 0.31415));
	real_vec2.store(std::vector<hmc_float>(5, 0.51413));
	for(int i=0; i<5; i++){
	  for(int j=0; j<5; j++){
	    physics::lattices::saxpby(&sf, real_vec1, i, cold, real_vec2, j, zero);
	    BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.04934511125, 1.e-8);
	    physics::lattices::saxpby(&sf, real_vec1, i, cold, real_vec2, j, cold);
	    BOOST_CHECK_CLOSE(physics::lattices::squarenorm(sf), 0.3430238792, 1.e-8);
	  }
	}
}

BOOST_AUTO_TEST_CASE(saxpbypz)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	Staggeredfield_eo gaussian(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	gaussian.set_gaussian(prng);
	Staggeredfield_eo cold(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	cold.set_cold();
	Staggeredfield_eo zero(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	zero.set_zero();
	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

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

BOOST_AUTO_TEST_CASE(sax_vec_and_squarenorm)
{
	using physics::lattices::Staggeredfield_eo;

	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);

	physics::lattices::Vector<hmc_float> zeros(3, system);
	physics::lattices::Vector<hmc_float> ones(3, system);
	physics::lattices::Vector<hmc_float> alpha(3, system);
	physics::lattices::Vector<hmc_float> result(3, system);
	std::vector<hmc_float> zeros_host(3, 0.);
	std::vector<hmc_float> ones_host(3, 1.);
	std::vector<hmc_float> alpha_host = {1., 1.25, 1.5};
	std::vector<hmc_float> reference = {0.5, 0.78125, 1.125};
	zeros.store(zeros_host);
	ones.store(ones_host);
	alpha.store(alpha_host);
	
	Staggeredfield_eo cold(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	cold.set_cold();
	Staggeredfield_eo rnd(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	physics::lattices::pseudo_randomize<Staggeredfield_eo, su3vec>(&rnd, 123);

	physics::lattices::sax_vec_and_squarenorm(&result, zeros, rnd);
	std::vector<hmc_float> result_host = result.get();
	for(uint i=0; i<result_host.size(); i++)
	  BOOST_CHECK_EQUAL(result_host[i], 0);
	
 	physics::lattices::sax_vec_and_squarenorm(&result, ones, cold);
 	result_host = result.get();
 	for(uint i=0; i<result_host.size(); i++)
 	  BOOST_CHECK_CLOSE(result_host[i], 0.5, 1.e-8);

	physics::lattices::sax_vec_and_squarenorm(&result, alpha, cold);
	result_host = result.get();
	for(uint i=0; i<result_host.size(); i++)
	  BOOST_CHECK_CLOSE(result_host[i], reference[i], 1.e-8);	
}

BOOST_AUTO_TEST_CASE(pseudorandomize)
{
	using physics::lattices::Staggeredfield_eo;
	
	const char * _params[] = {"foo", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters(params);
	physics::PRNG prng(system, &prngParameters);
	
	Staggeredfield_eo sf(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	sf.set_zero();
	BOOST_CHECK_EQUAL(physics::lattices::squarenorm(sf), 0);
	physics::lattices::pseudo_randomize<Staggeredfield_eo, su3vec>(&sf, 123);
	logger.info() << "The squarenorm of the pseudorandomized field is " << physics::lattices::squarenorm(sf);
}

