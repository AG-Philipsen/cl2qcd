/** @file
 * Tests of the explicit fermionmatrix implementations
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schaefer, Stefano Lottini, Alessandro Sciarra
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

#include "fermionmatrix.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::fermionmatrix::explicit
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(M_wilson)
{
//void M_wilson(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 1);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 2);

		physics::fermionmatrix::M_wilson(&sf2, gf, sf1, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2610.3804893063798, 0.01);
		physics::fermionmatrix::M_wilson(&sf1, gf, sf2, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4356.332327032359, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 2);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 1);

		physics::fermionmatrix::M_wilson(&sf2, gf, sf1, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2655.7059719467552, 0.01);
		physics::fermionmatrix::M_wilson(&sf1, gf, sf2, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4487.2227568771214, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_plus)
{
//void M_tm_plus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 3);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 1);

		physics::fermionmatrix::M_tm_plus(&sf2, gf, sf1, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2573.5343424130424, 0.01);
		physics::fermionmatrix::M_tm_plus(&sf1, gf, sf2, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4301.7240841150415, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 4);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 2);

		physics::fermionmatrix::M_tm_plus(&sf2, gf, sf1, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2536.3274644928938, 0.01);
		physics::fermionmatrix::M_tm_plus(&sf1, gf, sf2, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4202.4585412814668, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_minus)
{
//void M_tm_minus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 5);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 3);

		physics::fermionmatrix::M_tm_minus(&sf2, gf, sf1, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2580.3455858426491, 0.01);
		physics::fermionmatrix::M_tm_minus(&sf1, gf, sf2, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4323.4437436445069, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 6);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 4);

		physics::fermionmatrix::M_tm_minus(&sf2, gf, sf1, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2648.2135270529998, 0.01);
		physics::fermionmatrix::M_tm_minus(&sf1, gf, sf2, params.get_kappa(), meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4457.8676622761632, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_inverse_sitediagonal)
{
//void M_tm_inverse_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};

		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 5);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::M_tm_inverse_sitediagonal(&sf2, sf1, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 4133.0112721801961, 0.01);
		physics::fermionmatrix::M_tm_inverse_sitediagonal(&sf1, sf2, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4133.0019729257556, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_sitediagnoal)
{
//void M_tm_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};

		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 7);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::M_tm_sitediagonal(&sf2, sf1, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 4142.9167782423428, 0.01);
		physics::fermionmatrix::M_tm_sitediagonal(&sf1, sf2, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4142.9260998050931, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_inverse_sitediagonal_minus)
{
//void M_tm_inverse_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};

		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 9);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::M_tm_inverse_sitediagonal_minus(&sf2, sf1, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 4095.9297157806914, 0.01);
		physics::fermionmatrix::M_tm_inverse_sitediagonal_minus(&sf1, sf2, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4095.9204999595668, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(M_tm_sitediagonal_minus)
{
//void M_tm_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--fermact=twistedmass"};
		meta::Inputparameters params(3, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};

		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 11);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::M_tm_sitediagonal_minus(&sf2, sf1, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf2), 4099.5308263836096, 0.01);
		physics::fermionmatrix::M_tm_sitediagonal_minus(&sf1, sf2, meta::get_mubar(params));
		BOOST_CHECK_CLOSE(squarenorm(sf1), 4099.5400503279689, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(dslash)
{
//void dslash(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, int evenodd, hmc_float kappa = ARG_DEF);

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 13);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::dslash(&sf2, gf, sf1, EVEN, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 3311.2698428285048, 0.01);
		physics::fermionmatrix::dslash(&sf1, gf, sf2, ODD, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf1), 3146.1039504225546, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());

		pseudo_randomize<Spinorfield, spinor>(&src, 28);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::fermionmatrix::dslash(&sf2, gf, sf1, EVEN, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf2), 251.84231257415126, 0.01);
		physics::fermionmatrix::dslash(&sf1, gf, sf2, ODD, params.get_kappa());
		BOOST_CHECK_CLOSE(squarenorm(sf1), 75.926255640020059, 0.01);
	}

}
