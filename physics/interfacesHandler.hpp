/** @file
 * interfacesHandler declaration
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

#include "algorithms/algorithmsInterfaces.hpp"
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
    	class M;
    	class Qplus;
    	class Qminus;
    	class QplusQminus;
    	class Aee;
    	class Aee_minus;
    	class Qplus_eo;
    	class Qminus_eo;
    	class QplusQminus_eo;
    	class D_KS_eo;
    	class MdagM_eo;
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
    struct InterfaceType<physics::fermionmatrix::M> {
            using value = physics::fermionmatrix::FermionmatrixParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Qplus> {
            using value = physics::fermionmatrix::FermionmatrixParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Qminus> {
            using value = physics::fermionmatrix::FermionmatrixParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::QplusQminus> {
            using value = physics::FermionParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Aee> {
            using value = physics::FermionEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Aee_minus> {
            using value = physics::FermionEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Qplus_eo> {
            using value = physics::FermionEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::Qminus_eo> {
            using value = physics::FermionEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::QplusQminus_eo> {
            using value = physics::FermionEoParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::D_KS_eo> {
            using value = physics::fermionmatrix::FermionmatrixStaggeredParametersInterface;
    };
    template<>
    struct InterfaceType<physics::fermionmatrix::MdagM_eo> {
            using value = physics::FermionStaggeredEoParametersInterface;
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

            //NOTE: Here the interfaces for the algorithms. Since the algorithms are functions and not objects,
            //      these getters cannot be included in the template above (getInterface).
            virtual const physics::algorithms::SolversParametersInterface& getSolversParametersInterface() = 0;
            virtual const physics::algorithms::MinMaxEigenvalueParametersInterface& getMinMaxEigenvalueParametersInterface() = 0;
            virtual const physics::algorithms::ForcesParametersInterface& getForcesParametersInterface() = 0;
            virtual const physics::algorithms::InversionParemetersInterface& getInversionParemetersInterface() = 0;
            virtual const physics::algorithms::IntegratorParametersInterface& getIntegratorParametersInterface() = 0;
            virtual const physics::algorithms::MolecularDynamicsInterface& getMolecularDynamicsInterface() = 0;
            virtual const physics::algorithms::MetropolisParametersInterface& getMetropolisParametersInterface() = 0;
            virtual const physics::algorithms::HmcParametersInterface& getHmcParametersInterface() = 0;
            virtual const physics::algorithms::RhmcParametersInterface& getRhmcParametersInterface() = 0;

        private:
            virtual const physics::lattices::GaugefieldParametersInterface& getGaugefieldParametersInterface() = 0;
            virtual const physics::lattices::GaugemomentaParametersInterface& getGaugemomentaParametersInterface() = 0;
            virtual const physics::lattices::SpinorfieldParametersInterface& getSpinorfieldParametersInterface() = 0;
            virtual const physics::lattices::SpinorfieldEoParametersInterface& getSpinorfieldEoParametersInterface() = 0;
            virtual const physics::lattices::StaggeredfieldEoParametersInterface& getStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::lattices::RootedStaggeredfieldEoParametersInterface& getRootedStaggeredfieldEoParametersInterface() = 0;
            virtual const physics::fermionmatrix::FermionmatrixParametersInterface& getFermionmatrixParametersInterface() = 0;
            virtual const physics::fermionmatrix::FermionmatrixStaggeredParametersInterface& getFermionmatrixStaggeredParametersInterface() = 0;
            virtual const physics::FermionParametersInterface& getFermionParametersInterface() = 0;
            virtual const physics::FermionEoParametersInterface& getFermionEoParametersInterface() = 0;
            virtual const physics::FermionStaggeredEoParametersInterface& getFermionStaggeredEoParametersInterface() = 0;
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
    template<> inline const typename InterfaceType<physics::fermionmatrix::M>::value& InterfacesHandler::getInterface<physics::fermionmatrix::M>()
    {
        return getFermionmatrixParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Qplus>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Qplus>()
    {
        return getFermionmatrixParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Qminus>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Qminus>()
    {
        return getFermionmatrixParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::QplusQminus>::value& InterfacesHandler::getInterface<physics::fermionmatrix::QplusQminus>()
    {
        return getFermionParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Aee>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Aee>()
    {
        return getFermionEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Aee_minus>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Aee_minus>()
    {
        return getFermionEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Qplus_eo>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Qplus_eo>()
    {
        return getFermionEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::Qminus_eo>::value& InterfacesHandler::getInterface<physics::fermionmatrix::Qminus_eo>()
    {
        return getFermionEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::QplusQminus_eo>::value& InterfacesHandler::getInterface<physics::fermionmatrix::QplusQminus_eo>()
    {
        return getFermionEoParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::D_KS_eo>::value& InterfacesHandler::getInterface<physics::fermionmatrix::D_KS_eo>()
    {
        return getFermionmatrixStaggeredParametersInterface();
    }
    template<> inline const typename InterfaceType<physics::fermionmatrix::MdagM_eo>::value& InterfacesHandler::getInterface<physics::fermionmatrix::MdagM_eo>()
    {
        return getFermionStaggeredEoParametersInterface();
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
                  fermionmatrixStaggeredParametersInterface{nullptr},
                  fermionParametersInterface{nullptr},
                  fermionEoParametersInterface{nullptr},
                  fermionStaggeredEoParametersInterface{nullptr},
                  gaugeObservablesParametersInterface{nullptr},
                  wilsonTwoFlavourChiralCondensateParametersInterface{nullptr},
                  staggeredChiralCondensateParametersInterface{nullptr},
                  wilsonTwoFlavourCorrelatorsParametersInterface{nullptr},
                  solversParametersInterface{nullptr},
                  minMaxEigenvalueParametersInterface{nullptr},
                  forcesParametersInterface{nullptr},
                  inversionParametersInterface{nullptr},
                  integratorParametersInterface{nullptr},
                  molecularDynamicsParametersInterface{nullptr},
                  metropolisParametersInterface{nullptr},
                  hmcParametersInterface{nullptr},
                  rhmcParametersInterface{nullptr} {}
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
            const physics::algorithms::SolversParametersInterface& getSolversParametersInterface()
            {
                if( solversParametersInterface == nullptr)
                    solversParametersInterface = std::unique_ptr<const physics::algorithms::SolversParametersImplementation>(new physics::algorithms::SolversParametersImplementation{parameters});
                return *solversParametersInterface;
            }
            const physics::algorithms::MinMaxEigenvalueParametersInterface& getMinMaxEigenvalueParametersInterface()
            {
                if( minMaxEigenvalueParametersInterface == nullptr)
                    minMaxEigenvalueParametersInterface = std::unique_ptr<const physics::algorithms::MinMaxEigenvalueParametersImplementation>(new physics::algorithms::MinMaxEigenvalueParametersImplementation{parameters});
                return *minMaxEigenvalueParametersInterface;
            }
            const physics::algorithms::ForcesParametersInterface& getForcesParametersInterface()
            {
                if( forcesParametersInterface == nullptr)
                    forcesParametersInterface = std::unique_ptr<const physics::algorithms::ForcesParametersImplementation>(new physics::algorithms::ForcesParametersImplementation{parameters});
                return *forcesParametersInterface;
            }
            const physics::algorithms::InversionParemetersInterface& getInversionParemetersInterface()
            {
                if( inversionParametersInterface == nullptr)
                    inversionParametersInterface = std::unique_ptr<const physics::algorithms::InversionParametersImplementation>(new physics::algorithms::InversionParametersImplementation{parameters});
                return *inversionParametersInterface;
            }
            const physics::algorithms::IntegratorParametersInterface& getIntegratorParametersInterface()
            {
                if( integratorParametersInterface == nullptr)
                    integratorParametersInterface = std::unique_ptr<const physics::algorithms::IntegratorParametersImplementation>(new physics::algorithms::IntegratorParametersImplementation{parameters});
                return *integratorParametersInterface;
            }
            const physics::algorithms::MolecularDynamicsInterface& getMolecularDynamicsInterface()
            {
                if( molecularDynamicsParametersInterface == nullptr)
                    molecularDynamicsParametersInterface = std::unique_ptr<const physics::algorithms::MolecularDynamicsImplementation>(new physics::algorithms::MolecularDynamicsImplementation{parameters});
                return *molecularDynamicsParametersInterface;
            }
            const physics::algorithms::MetropolisParametersInterface& getMetropolisParametersInterface()
            {
                if( metropolisParametersInterface == nullptr)
                    metropolisParametersInterface = std::unique_ptr<const physics::algorithms::MetropolisParametersImplementation>(new physics::algorithms::MetropolisParametersImplementation{parameters});
                return *metropolisParametersInterface;
            }
            const physics::algorithms::HmcParametersInterface& getHmcParametersInterface()            {
                if( hmcParametersInterface == nullptr)
                     hmcParametersInterface = std::unique_ptr<const physics::algorithms::HmcParametersImplementation>(new physics::algorithms::HmcParametersImplementation{parameters});
                return *hmcParametersInterface;
            }
            const physics::algorithms::RhmcParametersInterface& getRhmcParametersInterface()
            {
                if( rhmcParametersInterface == nullptr)
                    rhmcParametersInterface = std::unique_ptr<const physics::algorithms::RhmcParametersImplementation>(new physics::algorithms::RhmcParametersImplementation{parameters});
                return *rhmcParametersInterface;
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
            const physics::fermionmatrix::FermionmatrixStaggeredParametersInterface& getFermionmatrixStaggeredParametersInterface() override
            {
                if(fermionmatrixStaggeredParametersInterface == nullptr)
                    fermionmatrixStaggeredParametersInterface = std::unique_ptr<const physics::fermionmatrix::FermionmatrixStaggeredParametersImplementation>(new physics::fermionmatrix::FermionmatrixStaggeredParametersImplementation{parameters});
                return *fermionmatrixStaggeredParametersInterface;
            }
            const physics::FermionParametersInterface& getFermionParametersInterface() override
            {
                if(fermionParametersInterface == nullptr)
                    fermionParametersInterface = std::unique_ptr<const physics::FermionParametersImplementation>(new physics::FermionParametersImplementation{parameters});
                return *fermionParametersInterface;
            }
            const physics::FermionEoParametersInterface& getFermionEoParametersInterface() override
            {
                if(fermionEoParametersInterface == nullptr)
                    fermionEoParametersInterface = std::unique_ptr<const physics::FermionEoParametersImplementation>(new physics::FermionEoParametersImplementation{parameters});
                return *fermionEoParametersInterface;
            }
            const physics::FermionStaggeredEoParametersInterface& getFermionStaggeredEoParametersInterface() override
            {
                if(fermionStaggeredEoParametersInterface == nullptr)
                    fermionStaggeredEoParametersInterface = std::unique_ptr<const physics::FermionStaggeredEoParametersImplementation>(new physics::FermionStaggeredEoParametersImplementation{parameters});
                return *fermionStaggeredEoParametersInterface;
            }

            const meta::Inputparameters& parameters;
            std::unique_ptr<const physics::lattices::GaugefieldParametersInterface> gaugefieldParametersInterface;
            std::unique_ptr<const physics::lattices::GaugemomentaParametersInterface> gaugemomentaParametersInterface;
            std::unique_ptr<const physics::lattices::SpinorfieldParametersInterface> spinorfieldParametersInterface;
            std::unique_ptr<const physics::lattices::SpinorfieldEoParametersInterface> spinorfieldEoParametersInterface;
            std::unique_ptr<const physics::lattices::StaggeredfieldEoParametersInterface> staggeredfieldEoParametersInterface;
            std::unique_ptr<const physics::lattices::RootedStaggeredfieldEoParametersInterface> rootedStaggaredfieldEoParametersInterface;
            std::unique_ptr<const physics::fermionmatrix::FermionmatrixParametersInterface> fermionmatrixParametersInterface;
            std::unique_ptr<const physics::fermionmatrix::FermionmatrixStaggeredParametersInterface> fermionmatrixStaggeredParametersInterface;
            std::unique_ptr<const physics::FermionParametersInterface> fermionParametersInterface;
            std::unique_ptr<const physics::FermionEoParametersInterface> fermionEoParametersInterface;
            std::unique_ptr<const physics::FermionStaggeredEoParametersInterface> fermionStaggeredEoParametersInterface;
            std::unique_ptr<const physics::observables::GaugeObservablesParametersInterface> gaugeObservablesParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface> wilsonTwoFlavourChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::StaggeredChiralCondensateParametersInterface> staggeredChiralCondensateParametersInterface;
            std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface> wilsonTwoFlavourCorrelatorsParametersInterface;
            std::unique_ptr<const physics::algorithms::SolversParametersInterface> solversParametersInterface;
            std::unique_ptr<const physics::algorithms::MinMaxEigenvalueParametersInterface> minMaxEigenvalueParametersInterface;
            std::unique_ptr<const physics::algorithms::ForcesParametersInterface> forcesParametersInterface;
            std::unique_ptr<const physics::algorithms::InversionParemetersInterface> inversionParametersInterface;
            std::unique_ptr<const physics::algorithms::IntegratorParametersInterface> integratorParametersInterface;
            std::unique_ptr<const physics::algorithms::MolecularDynamicsInterface> molecularDynamicsParametersInterface;
            std::unique_ptr<const physics::algorithms::MetropolisParametersInterface> metropolisParametersInterface;
            std::unique_ptr<const physics::algorithms::HmcParametersInterface> hmcParametersInterface;
            std::unique_ptr<const physics::algorithms::RhmcParametersInterface> rhmcParametersInterface;
    };

}
