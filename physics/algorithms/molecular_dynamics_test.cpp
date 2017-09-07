/** @file
 * Tests of the molecular dynamics algorithms
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

#include "molecular_dynamics.hpp"
#include "forces.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::molecular_dynamics
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../observables/gaugeObservables.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(md_update_gaugefield)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		{
			Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
			Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
			gm.zero();

			hmc_float ref = physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface());
			gauge_force(&gm, gf);
			physics::algorithms::md_update_gaugefield(&gf, gm, .5);
			double plaq = physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface());
			BOOST_CHECK_CLOSE(plaq, ref, 0.01);
		}

		{
			Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
			Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
			pseudo_randomize<Gaugemomenta, ae>(&gm, 415);

			physics::algorithms::md_update_gaugefield(&gf, gm, .5);
			double plaq = physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface());
			BOOST_CHECK_CLOSE(plaq, 0.80918156710730049, 0.01);
		}
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		pseudo_randomize<Gaugemomenta, ae>(&gm, 123);

		double plaq = physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface());
		BOOST_CHECK_CLOSE(plaq, 0.57107711169452713, 0.0001);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		plaq = physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface());
		BOOST_REQUIRE_CLOSE(plaq, 0.32089465123266286, 0.01);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		BOOST_REQUIRE_CLOSE(physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface()), 0.57107711169452713, 0.0001);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		BOOST_REQUIRE_CLOSE(physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface()), 0.57107711169452713, 0.0001);
		gauge_force(&gm, gf);

		BOOST_REQUIRE_CLOSE(squarenorm(gm), 52723.299867438494, 0.0001);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		BOOST_REQUIRE_CLOSE(physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface()), 0.0060440132434446334, 0.01);
		gauge_force(&gm, gf);
		BOOST_REQUIRE_CLOSE(squarenorm(gm), 82900.801546488685, 0.01);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		BOOST_REQUIRE_CLOSE(physics::observables::measurePlaquette(&gf, interfacesHandler.getGaugeObservablesParametersInterface()), -0.0076685322051177783, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 11);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 12);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2523.6197747176057, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 13);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 14);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2588.2852881545778, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 15);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1122.194230730885, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 16);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1114.3019247079062, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(mdUpdateStaggeredRootedSpinorfieldEo)
{
	using namespace physics::lattices;
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
	meta::Inputparameters params(4, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	bool inv = 0;
	physics::algorithms::Rational_Approximation approx(10,1,4,1e-5,1,inv);
	Rooted_Staggeredfield_eo in(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>());
	Rooted_Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);

	pseudo_randomize<Staggeredfield_eo, su3vec>(in[0].get(), 13);
	physics::algorithms::md_update_spinorfield(&out, gf, in, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());

	BOOST_CHECK_CLOSE(squarenorm(*out[0].get()), 328.75052967496629, 1e-8);
}

BOOST_AUTO_TEST_CASE(mdUpdateStaggeredRootedSpinorfieldEoWithPseudofermions)
{
	using namespace physics::lattices;
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1", "--num_pseudofermions=2"};
	meta::Inputparameters params(5, _params);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	hardware::HardwareParametersImplementation hP(&params);
	hardware::code::OpenClKernelParametersImplementation kP(params);
	hardware::System system(hP, kP);
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	bool inv = 0;
	physics::algorithms::Rational_Approximation approx(10,1,4,1e-5,1,inv);
	Rooted_Staggeredfield_eo in(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>());
	Rooted_Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);

	{
        for(const auto& in_j : in)
            pseudo_randomize<Staggeredfield_eo, su3vec>(in_j.get(), 13);
        physics::algorithms::md_update_spinorfield(&out, gf, in, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());

        for(const auto& out_j : out)
            BOOST_CHECK_CLOSE(squarenorm(*out_j.get()), 328.75052967496629, 1e-8);
    }
    {
        for(const auto& in_j : in)
            in_j.get()->set_cold();
        physics::algorithms::md_update_spinorfield(&out, gf, in, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());

        for(const auto& out_j : out)
            BOOST_CHECK_CLOSE(squarenorm(*out_j.get()), 0.66024054611635885, 1e-8);
    }
    {
        for(const auto& in_j : in)
            in_j.get()->set_gaussian(prng);
        physics::algorithms::md_update_spinorfield(&out, gf, in, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());

        BOOST_CHECK_CLOSE(squarenorm(*out[0].get()), 455.28018993352305, 1e-8);
        BOOST_CHECK_CLOSE(squarenorm(*out[1].get()), 538.11156457534264, 1e-8);
    }
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_mp)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 22);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2539.2078579177487, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 23);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 24);
		gm.zero();

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2658.5475465525888, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_mp_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 25);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1105.9891302664669, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 26);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1081.4370785519473, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_gauge)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_gauge(&gm, .5, gf, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_gauge(&gm, .5, gf, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 13180.824966859615, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_fermion)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 33);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_fermion(&gm, .6, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 16956.313363328729, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_fermion_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 25);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_fermion(&gm, .6, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 1324.8991304359643, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_detratio)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--kappa_mp=.25"};
		meta::Inputparameters params(3, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 35);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_detratio(&gm, .4, gf, sf1, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 2070.2139029781247, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_detratio_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--kappa_mp=.25"};
		meta::Inputparameters params(3, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 27);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_detratio(&gm, .4, gf, sf1, system, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 5835.9232744740457, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 37);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum(&gm, .5, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 22735.531268100596, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 29);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum(&gm, .6, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 20459.049258021565, 0.01);
	}
}
