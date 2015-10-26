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
#include "../executables/exceptions.h"

namespace physics {

    namespace lattices {
        class Gaugefield;
        class Gaugemomenta;
        class Spinorfield;
        class Spinorfield_eo;
        class Staggeredfield_eo;
        class Rooted_Staggeredfield_eo;
    }

    namespace fermionmatrix {
    	class Fermionmatrix;
    }

    template<typename NOT_IMPORTANT> struct InterfaceType; //Not defined in order to prevent non specialized template instatiation!

    template<>
    struct InterfaceType<physics::lattices::Gaugefield> {
            using value = physics::lattices::GaugefieldParametersInterface;
    };
    template<>
    struct InterfaceType<physics::lattices::Gaugemomenta> {
            using value = physics::lattices::GaugemomentaParametersInterface;
    };
    template<>
    struct InterfaceType<physics::lattices::Spinorfield> {
            using value = physics::lattices::SpinorfieldParametersInterface;
    };
    template<>
    struct InterfaceType<physics::lattices::Spinorfield_eo> {
            using value = physics::lattices::SpinorfieldEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::lattices::Staggeredfield_eo> {
            using value = physics::lattices::StaggeredfieldEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::lattices::Rooted_Staggeredfield_eo> {
            using value = physics::lattices::RootedStaggeredfieldEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Fermionmatrix> {
            using value = physics::fermionmatrix::FermionmatrixParametersInterface;
    };

    class InterfacesHandler {
        public:
            virtual ~InterfacesHandler() {};
            template<class OBJECT> const typename InterfaceType<OBJECT>::value& getInterface();//Not defined in order to prevent non specialized template instatiation!

            //TODO: For the moment there is no object responsible to calculate observables, so we leave the following getters out of the above template
            //      It would be better to create objects in physics::observables and move these getters in the template above.
            virtual const physics::observables::GaugeObservablesParametersInterface& getGaugeObservablesParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& getWilsonTwoFlavourChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::StaggeredChiralCondensateParametersInterface& getStaggeredChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsCondensateParametersInterface() = 0;
        private:
            virtual const physics::lattices::GaugefieldParametersInterface& getGaugefieldParametersInterface() = 0;
            virtual const physics::lattices::GaugemomentaParametersInterface& getGaugemomentaParametersInterface() = 0;
            virtual const physics::lattices::SpinorfieldParametersInterface& getSpinorfieldParametersInterface() = 0;
            virtual const physics::lattices::SpinorfieldEoParametersInterface& getSpinorfieldEoParametersInterface() = 0;
            virtual const physics::lattices::StaggeredfieldEoParametersInterface& getStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::lattices::RootedStaggeredfieldEoParametersInterface& getRootedStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::fermionmatrix::FermionmatrixParametersInterface& getFermionmatrixParametersInterface() = 0;
    };

    template<> inline const typename InterfaceType<physics::lattices::Gaugefield>::value& InterfacesHandler::getInterface<physics::lattices::Gaugefield>()
    {
        return getGaugefieldParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::lattices::Gaugemomenta>::value& InterfacesHandler::getInterface<physics::lattices::Gaugemomenta>()
    {
        return getGaugemomentaParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::lattices::Spinorfield>::value& InterfacesHandler::getInterface<physics::lattices::Spinorfield>()
    {
        return getSpinorfieldParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::lattices::Spinorfield_eo>::value& InterfacesHandler::getInterface<physics::lattices::Spinorfield_eo>()
    {
        return getSpinorfieldEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::lattices::Staggeredfield_eo>::value& InterfacesHandler::getInterface<physics::lattices::Staggeredfield_eo>()
    {
        return getStaggeredfieldEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::lattices::Rooted_Staggeredfield_eo>::value& InterfacesHandler::getInterface<physics::lattices::Rooted_Staggeredfield_eo>()
    {
        return getRootedStaggeredfieldEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Fermionmatrix>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Fermionmatrix>()
    {
        return getFermionmatrixParametersInterface();
    }

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
				  spinorfieldEoParametersInterface{nullptr},
                  staggeredfieldEoParametersInterface{nullptr},
                  rootedStaggaredfieldEoParametersInterface{nullptr},
                  fermionmatrixParametersInterface{nullptr},
                  gaugeObservablesParametersInterface{nullptr},
                  wilsonTwoFlavourChiralCondensateParametersInterface{nullptr},
                  staggeredChiralCondensateParametersInterface{nullptr},
                  wilsonTwoFlavourCorrelatorsParametersInterface{nullptr} {}
            ~InterfacesHandlerImplementation() {}
            const physics::observables::GaugeObservablesParametersInterface& getGaugeObservablesParametersInterface() override
            {
                if(gaugeObservablesParametersInterface == nullptr)
                    gaugeObservablesParametersInterface = std::unique_ptr<const physics::observables::GaugeObservablesParametersImplementation>(new physics::observables::GaugeObservablesParametersImplementation{parameters});
                return *gaugeObservablesParametersInterface;
            }
            const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& getWilsonTwoFlavourChiralCondensateParametersInterface() override
            {
                if(wilsonTwoFlavourChiralCondensateParametersInterface == nullptr)
                    wilsonTwoFlavourChiralCondensateParametersInterface = std::unique_ptr<const physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation>(new physics::observables::WilsonTwoFlavourChiralCondensateParametersImplementation {parameters});
                return *wilsonTwoFlavourChiralCondensateParametersInterface;
            }
            const physics::observables::StaggeredChiralCondensateParametersInterface& getStaggeredChiralCondensateParametersInterface() override
            {
                if(staggeredChiralCondensateParametersInterface == nullptr)
                    staggeredChiralCondensateParametersInterface = std::unique_ptr<const physics::observables::StaggeredChiralCondensateParametersImplementation>(new physics::observables::StaggeredChiralCondensateParametersImplementation{parameters});
                return *staggeredChiralCondensateParametersInterface;
            }
            const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsCondensateParametersInterface() override
            {
                if(wilsonTwoFlavourCorrelatorsParametersInterface == nullptr)
                    wilsonTwoFlavourCorrelatorsParametersInterface = std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation>(new physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation{parameters});
                return *wilsonTwoFlavourCorrelatorsParametersInterface;
            }

        private:
            const physics::lattices::GaugefieldParametersInterface& getGaugefieldParametersInterface() override
            {
                if(gaugefieldParametersInterface == nullptr)
                    gaugefieldParametersInterface = std::unique_ptr<const physics::lattices::GaugefieldParametersImplementation>(new physics::lattices::GaugefieldParametersImplementation{&parameters});
                return *gaugefieldParametersInterface;
            }
            const physics::lattices::GaugemomentaParametersInterface& getGaugemomentaParametersInterface() override
            {
                if(gaugemomentaParametersInterface == nullptr)
                    gaugemomentaParametersInterface = std::unique_ptr<const physics::lattices::GaugemomentaParametersImplementation>(new physics::lattices::GaugemomentaParametersImplementation{parameters});
                return *gaugemomentaParametersInterface;
            }
            const physics::lattices::SpinorfieldParametersInterface& getSpinorfieldParametersInterface() override
            {
                if(spinorfieldParametersInterface == nullptr)
                    spinorfieldParametersInterface = std::unique_ptr<const physics::lattices::SpinorfieldParametersImplementation>(new physics::lattices::SpinorfieldParametersImplementation{parameters});
                return *spinorfieldParametersInterface;
            }
            const physics::lattices::SpinorfieldEoParametersInterface& getSpinorfieldEoParametersInterface() override
            {
            	if(spinorfieldEoParametersInterface == nullptr)
            		spinorfieldEoParametersInterface = std::unique_ptr<const physics::lattices::SpinorfieldEoParametersImplementation>(new physics::lattices::SpinorfieldEoParametersImplementation{parameters});
            	return *spinorfieldEoParametersInterface;
            }
            const physics::lattices::StaggeredfieldEoParametersInterface& getStaggeredfieldEoParametersInterface() override
            {
                if(staggeredfieldEoParametersInterface == nullptr)
                    staggeredfieldEoParametersInterface = std::unique_ptr<const physics::lattices::StaggeredfieldEoParametersImplementation>(new physics::lattices::StaggeredfieldEoParametersImplementation{parameters});
                return *staggeredfieldEoParametersInterface;
            }
            const physics::lattices::RootedStaggeredfieldEoParametersInterface& getRootedStaggeredfieldEoParametersInterface() override
            {
                if(rootedStaggaredfieldEoParametersInterface == nullptr)
                    rootedStaggaredfieldEoParametersInterface = std::unique_ptr<const physics::lattices::RootedStaggeredfieldEoParametersImplementation>(new physics::lattices::RootedStaggeredfieldEoParametersImplementation{parameters});
                return *rootedStaggaredfieldEoParametersInterface;
            }
            const physics::fermionmatrix::FermionmatrixParametersInterface& getFermionmatrixParametersInterface() override
            {
                if(fermionmatrixParametersInterface == nullptr)
                    fermionmatrixParametersInterface = std::unique_ptr<const physics::fermionmatrix::FermionmatrixParametersImplementation>(new physics::fermionmatrix::FermionmatrixParametersImplementation{parameters});
                return *fermionmatrixParametersInterface;
            }

            const meta::Inputparameters& parameters;
            std::unique_ptr<const physics::lattices::GaugefieldParametersInterface> gaugefieldParametersInterface;
            std::unique_ptr<const physics::lattices::GaugemomentaParametersInterface> gaugemomentaParametersInterface;
            std::unique_ptr<const physics::lattices::SpinorfieldParametersInterface> spinorfieldParametersInterface;
            std::unique_ptr<const physics::lattices::SpinorfieldEoParametersInterface> spinorfieldEoParametersInterface;
            std::unique_ptr<const physics::lattices::StaggeredfieldEoParametersInterface> staggeredfieldEoParametersInterface;
            std::unique_ptr<const physics::lattices::RootedStaggeredfieldEoParametersInterface> rootedStaggaredfieldEoParametersInterface;
            std::unique_ptr<const physics::fermionmatrix::FermionmatrixParametersInterface> fermionmatrixParametersInterface;
            std::unique_ptr<const physics::observables::GaugeObservablesParametersInterface> gaugeObservablesParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface> wilsonTwoFlavourChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::StaggeredChiralCondensateParametersInterface> staggeredChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface> wilsonTwoFlavourCorrelatorsParametersInterface;
    };

}
