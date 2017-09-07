/** @file
 * Tests of the metropolis algorithms
 *
 * Copyright (c) 2017 Francesca Cuteri <cuteri@th.physik.uni-frankfurt.de>
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

#include "metropolis.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::molecular_dynamics
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(metropolisStaggeredRootedSpinorfieldEo)
{
    using namespace physics::lattices;
    physics::algorithms::Rational_Approximation approx(8, 1,2, 1.e-5,1);
    const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
    meta::Inputparameters params(4, _params);
    physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
    physics::PrngParametersImplementation prngParameters{params};
    physics::PRNG prng{system, &prngParameters};

    Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);

    Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

    Rooted_Staggeredfield_eo phi(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);

    hmc_float beta = 6.0;
    hmc_float spinor_energy_init = 0.;
    gm.zero();

    {
        phi[0].get()->set_zero();
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_EQUAL(0,obs.deltaH);
    }

    {
        phi[0].get()->set_cold();
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_CLOSE(obs.deltaH, -4.999853415117, 0.01);
    }

    {
        phi[0].get()->set_gaussian(prng);
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_CLOSE(obs.deltaH, -570.2704806215952, 0.01);
    }
}

BOOST_AUTO_TEST_CASE(metropolisStaggeredRootedSpinorfieldEoWithPseudofermions)
{
    using namespace physics::lattices;
    physics::algorithms::Rational_Approximation approx(8, 1,2, 1.e-5,1);
    const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1", "--num_pseudofermions=2"};
    meta::Inputparameters params(5, _params);
    physics::InterfacesHandlerImplementation interfacesHandler{params};
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
    physics::PrngParametersImplementation prngParameters{params};
    physics::PRNG prng{system, &prngParameters};

    Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);

    Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

    Rooted_Staggeredfield_eo phi(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);

    hmc_float beta = 6.0;
    hmc_float spinor_energy_init = 0.;
    gm.zero();

    {
        for(const auto& phi_j : phi)
            phi_j.get()->set_zero();
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_EQUAL(0,obs.deltaH);
    }

    {
        for(const auto& phi_j : phi)
            phi_j.get()->set_cold();
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_CLOSE(obs.deltaH, -4.999853415117*2, 0.01);
        //the reference value here is double the one of the metropolisStaggeredRootedSpinorfieldEo TEST_CASE with the same cold initialization of the spinorfield since we use the same Rational_Approximation
    }

    {
        for(const auto& phi_j : phi)
            phi_j.get()->set_gaussian(prng);
        const hmc_observables obs = physics::algorithms::metropolis(0., beta, gf, gf, gm, gm, phi, spinor_energy_init,
                                                              nullptr, spinor_energy_init, system, interfacesHandler);
        BOOST_CHECK_CLOSE(obs.deltaH, -1059.318832641, 0.01);
    }
}
