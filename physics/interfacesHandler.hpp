/** @file
 * PRNG PRNG unit declaration
 *
 * Copyright 2015 Alessandro Sciarra, Christopher Czaban
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

#pragma once

#include "fermionmatrix/fermionmatrixInterfaces.hpp"
#include "lattices/latticesInterfaces.hpp"
#include "observables/observablesInterfaces.hpp"

namespace physics {

    class InterfacesHandler {
        public:
            virtual ~InterfacesHandler() {};
            //The following getters are not constant since in the implementation they initialize a private member the first time they are used!
            virtual const physics::lattices::GaugefieldParametersInterface& getGaugefieldParametersInterface() = 0;
            virtual const physics::lattices::GaugemomentaParametersInterface& getGaugemomentaParametersInterface() = 0;
            virtual const physics::lattices::SpinorfieldParametersInterface& getSpinorfieldParametersInterface() = 0;
            virtual const physics::lattices::StaggaredfieldEoParametersInterface& getStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::lattices::RootedStaggaredfieldEoParametersInterface& getRootedStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::fermionmatrix::FermionmatrixParametersInterface& getFermionmatrixParametersInterface() = 0;
            virtual const physics::observables::GaugeObservablesParametersInterface& getGaugeObservablesParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& getWilsonTwoFlavourChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::StaggeredChiralCondensateParametersInterface& getStaggeredChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsCondensateParametersInterface() = 0;
    };

}

#include <memory>

namespace physics {

    class InterfacesHandlerImplementation final : public InterfacesHandler {
        public:
            InterfacesHandlerImplementation() = delete;
            InterfacesHandlerImplementation(const meta::Inputparameters& parametersIn)
                : parameters(parametersIn),
                  gaugefieldParametersInterface{nullptr},
                  gaugemomentaParametersInterface{nullptr},
                  spinorfieldParametersInterface{nullptr},
                  staggaredfieldEoParametersInterface{nullptr},
                  rootedStaggaredfieldEoParametersInterface{nullptr},
                  fermionmatrixParametersInterface{nullptr},
                  gaugeObservablesParametersInterface{nullptr},
                  wilsonTwoFlavourChiralCondensateParametersInterface{nullptr},
                  staggeredChiralCondensateParametersInterface{nullptr},
                  wilsonTwoFlavourCorrelatorsParametersInterface{nullptr} {}
            ~InterfacesHandlerImplementation() {}
            const physics::lattices::GaugefieldParametersInterface& getGaugefieldParametersInterface()
            {
                if(gaugefieldParametersInterface == nullptr)
                    gaugefieldParametersInterface = std::unique_ptr<const physics::lattices::GaugefieldParametersImplementation>(new physics::lattices::GaugefieldParametersImplementation{&parameters});
                return *gaugefieldParametersInterface;
            }
            const physics::lattices::GaugemomentaParametersInterface& getGaugemomentaParametersInterface()
            {
                if(gaugemomentaParametersInterface == nullptr)
                    gaugemomentaParametersInterface = std::unique_ptr<const physics::lattices::GaugemomentaParametersImplementation>(new physics::lattices::GaugemomentaParametersImplementation{parameters});
                return *gaugemomentaParametersInterface;
            }
            const physics::lattices::SpinorfieldParametersInterface& getSpinorfieldParametersInterface()
            {
                if(spinorfieldParametersInterface == nullptr)
                    spinorfieldParametersInterface = std::unique_ptr<const physics::lattices::SpinorfieldParametersImplementation>(new physics::lattices::SpinorfieldParametersImplementation{parameters});
                return *spinorfieldParametersInterface;
            }
            const physics::lattices::StaggaredfieldEoParametersInterface& getStaggeredfieldEoParametersInterface()
            {
                if(staggaredfieldEoParametersInterface == nullptr)
                    staggaredfieldEoParametersInterface = std::unique_ptr<const physics::lattices::StaggaredfieldEoParametersImplementation>(new physics::lattices::StaggaredfieldEoParametersImplementation{parameters});
                return *staggaredfieldEoParametersInterface;
            }
            const physics::lattices::RootedStaggaredfieldEoParametersInterface& getRootedStaggeredfieldEoParametersInterface()
            {
                if(rootedStaggaredfieldEoParametersInterface == nullptr)
                    rootedStaggaredfieldEoParametersInterface = std::unique_ptr<const physics::lattices::RootedStaggaredfieldEoParametersImplementation>(new physics::lattices::RootedStaggaredfieldEoParametersImplementation{parameters});
                return *rootedStaggaredfieldEoParametersInterface;
            }
            const physics::fermionmatrix::FermionmatrixParametersInterface& getFermionmatrixParametersInterface()
            {
                if(fermionmatrixParametersInterface == nullptr)
                    fermionmatrixParametersInterface = std::unique_ptr<const physics::fermionmatrix::FermionmatrixParametersImplementation>(new physics::fermionmatrix::FermionmatrixParametersImplementation{parameters});
                return *fermionmatrixParametersInterface;
            }
            const physics::observables::GaugeObservablesParametersInterface& getGaugeObservablesParametersInterface()
            {
                if(gaugeObservablesParametersInterface == nullptr)
                    gaugeObservablesParametersInterface = std::unique_ptr<const physics::observables::GaugeObservablesParametersImplementation>(new physics::observables::GaugeObservablesParametersImplementation{parameters});
                return *gaugeObservablesParametersInterface;
            }
            const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& getWilsonTwoFlavourChiralCondensateParametersInterface()
            {
                if(wilsonTwoFlavourChiralCondensateParametersInterface == nullptr)
                    wilsonTwoFlavourChiralCondensateParametersInterface = std::unique_ptr<const physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation>(new physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation {parameters});
                return *wilsonTwoFlavourChiralCondensateParametersInterface;
            }
            const physics::observables::StaggeredChiralCondensateParametersInterface& getStaggeredChiralCondensateParametersInterface()
            {
                if(staggeredChiralCondensateParametersInterface == nullptr)
                    staggeredChiralCondensateParametersInterface = std::unique_ptr<const physics::observables::StaggeredChiralCondensateParametersImplementation>(new physics::observables::StaggeredChiralCondensateParametersImplementation{parameters});
                return *staggeredChiralCondensateParametersInterface;
            }
            const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsCondensateParametersInterface()
            {
                if(wilsonTwoFlavourCorrelatorsParametersInterface == nullptr)
                    wilsonTwoFlavourCorrelatorsParametersInterface = std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation>(new physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation{parameters});
                return *wilsonTwoFlavourCorrelatorsParametersInterface;
            }

        private:
            const meta::Inputparameters& parameters;
            std::unique_ptr<const physics::lattices::GaugefieldParametersInterface> gaugefieldParametersInterface;
            std::unique_ptr<const physics::lattices::GaugemomentaParametersInterface> gaugemomentaParametersInterface;
            std::unique_ptr<const physics::lattices::SpinorfieldParametersInterface> spinorfieldParametersInterface;
            std::unique_ptr<const physics::lattices::StaggaredfieldEoParametersInterface> staggaredfieldEoParametersInterface;
            std::unique_ptr<const physics::lattices::RootedStaggaredfieldEoParametersInterface> rootedStaggaredfieldEoParametersInterface;
            std::unique_ptr<const physics::fermionmatrix::FermionmatrixParametersInterface> fermionmatrixParametersInterface;
            std::unique_ptr<const physics::observables::GaugeObservablesParametersInterface> gaugeObservablesParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface> wilsonTwoFlavourChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::StaggeredChiralCondensateParametersInterface> staggeredChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface> wilsonTwoFlavourCorrelatorsParametersInterface;
    };

}
