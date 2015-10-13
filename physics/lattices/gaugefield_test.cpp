/** @file
 * Unit test for the physics::lattices::Gaugefield class
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

#include "gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugefield
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>

#include "../observables/gaugeObservables.h"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	{
		const char * _params[] = {"foo", "--ntime=8"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		// init hot
		Gaugefield gf2(system, &parametersTmp, prng, true);

		// init cold
		Gaugefield gf3(system, &parametersTmp, prng, false);
		logger.fatal() << physics::observables::measurePlaquette(&gf3);
		BOOST_CHECK_CLOSE(physics::observables::measurePlaquette(&gf3), 1., 0.1);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		// init from file
		Gaugefield gf(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/hardware/code/conf.00200");
		BOOST_CHECK_CLOSE(physics::observables::measurePlaquette(&gf), 0.57107711169452713, 0.1);
	}
}

void test_save(bool hot) {
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	const GaugefieldParametersImplementation parametersTmp{ &params };
	hardware::System system(params);
	physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
	physics::PRNG prng(system, &prngParameters);

	Gaugefield gf(system, &parametersTmp, prng, hot);
	gf.save("conf.test", 0);

	hmc_float orig_plaq, reread_plaq;
	hmc_complex orig_pol, reread_pol;

	orig_plaq = physics::observables::measurePlaquette(&gf);
	orig_pol = physics::observables::measurePolyakovloop(&gf);

	//NOTE: the conversion to std::string is necessary, otherwise the compiler creates a boolean!
	Gaugefield reread(system, &parametersTmp, prng, (std::string) "conf.test");
	
	reread_plaq =  physics::observables::measurePlaquette(&reread);
	reread_pol = physics::observables::measurePolyakovloop(&reread);

	BOOST_CHECK_EQUAL(orig_plaq, reread_plaq);
	BOOST_CHECK_EQUAL(orig_pol, reread_pol);
}

BOOST_AUTO_TEST_CASE(save)
{
	test_save(false);
	test_save(true);
}

BOOST_AUTO_TEST_CASE(rectangles)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--gaugeact=wilson", "--ntime=4"};
	meta::Inputparameters params(3, _params);
	const GaugefieldParametersImplementation parametersTmp{ &params };
	hardware::System system(params);
	physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
	physics::PRNG prng(system, &prngParameters);

	Gaugefield gf(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/hardware/code/conf.00200");
	BOOST_CHECK_THROW( physics::observables::measureRectangles(&gf);, std::logic_error);

	const char * _params2[] = {"foo", "--gaugeact=tlsym", "--ntime=4"};
	meta::Inputparameters params2(3, _params2);
	const GaugefieldParametersImplementation parametersTmp2{ &params2 };
	hardware::System system2(params2);
	physics::ParametersPrng_fromMetaInputparameters prngParameters2(&params2);
	physics::PRNG prng2(system, &prngParameters2);

	Gaugefield gf2(system2, &parametersTmp2, prng2, std::string(SOURCEDIR) + "/hardware/code/conf.00200");
	BOOST_CHECK_CLOSE(physics::observables::measureRectangles(&gf2), 1103.2398401620451, 0.1);
}

BOOST_AUTO_TEST_CASE(polyakov)
{
	using namespace physics::lattices;

	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		Gaugefield gf(system, &parametersTmp, prng, false);

		hmc_complex pol = physics::observables::measurePolyakovloop(&gf);;
		BOOST_CHECK_CLOSE(pol.re, 1., 0.1);
		BOOST_CHECK_CLOSE(pol.im, 0., 0.1);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		Gaugefield gf(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/hardware/code/conf.00200");
		hmc_complex pol = physics::observables::measurePolyakovloop(&gf);
		BOOST_CHECK_CLOSE(pol.re, -0.11349672123636857, 0.1);
		BOOST_CHECK_CLOSE(pol.im, 0.22828243566855227, 0.1);
	}
}

BOOST_AUTO_TEST_CASE(halo_update)
{
	using namespace physics::lattices;

	hmc_float orig_plaq, new_plaq;
	hmc_complex orig_pol, new_pol;

	// simple test, gaugeobservables should not get changed by halo exchange
	// if the original gaugefield is given
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		Gaugefield gf(system, &parametersTmp, prng, false);

		orig_plaq = physics::observables::measurePlaquette(&gf);
		orig_pol = physics::observables::measurePolyakovloop(&gf);

		gf.update_halo();

		new_plaq = physics::observables::measurePlaquette(&gf);
		new_pol = physics::observables::measurePolyakovloop(&gf);

		BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
		BOOST_CHECK_EQUAL(orig_pol, new_pol);
	}

	{
		const char * _params[] = {"foo"};
		meta::Inputparameters params(1, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		Gaugefield gf(system, &parametersTmp, prng, true);

		orig_plaq = physics::observables::measurePlaquette(&gf);
		orig_pol = physics::observables::measurePolyakovloop(&gf);

		gf.update_halo();

		new_plaq = physics::observables::measurePlaquette(&gf);
		new_pol = physics::observables::measurePolyakovloop(&gf);

		BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
		BOOST_CHECK_EQUAL(orig_pol, new_pol);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		const GaugefieldParametersImplementation parametersTmp{ &params };
		hardware::System system(params);
		physics::ParametersPrng_fromMetaInputparameters prngParameters(&params);
		physics::PRNG prng(system, &prngParameters);

		Gaugefield gf(system, &parametersTmp, prng, std::string(SOURCEDIR) + "/hardware/code/conf.00200");

		orig_plaq = physics::observables::measurePlaquette(&gf);
		orig_pol = physics::observables::measurePolyakovloop(&gf);

		gf.update_halo();

		new_plaq = physics::observables::measurePlaquette(&gf);
		new_pol = physics::observables::measurePolyakovloop(&gf);

		BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
		BOOST_CHECK_EQUAL(orig_pol, new_pol);
	}
}
