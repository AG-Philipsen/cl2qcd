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

#include "additionalParameters.hpp"
#include "sourcesInterface.hpp"
#include "algorithms/algorithmsInterfaces.hpp"
#include "fermionmatrix/fermionmatrixInterfaces.hpp"
#include "lattices/latticesInterfaces.hpp"
#include "observables/observablesInterfaces.hpp"

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

    //The following templates are not defined in order to prevent non specialized template instantiation!
    template<typename NOT_IMPORTANT> struct InterfaceType;

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
            //The following templates are not defined in order to prevent non specialized template instatiation!
            template<class OBJECT> const typename InterfaceType<OBJECT>::value& getInterface();
            template<class OBJECT> const physics::AdditionalParameters& getAdditionalParameters(bool withMassPreconditioning = false);

            //TODO: For the moment there is no object responsible to calculate observables, so we leave the following getters out of the above template
            //      Think about to create objects in physics::observables and move these getters in the template above.
            //      Actually the observables are similar to algorithms in the sense that their nature is more a function than an object.
            //      Nevertheless, one could make them static in a class so that one would have here an object to use in getInterface!
            virtual const physics::observables::GaugeObservablesParametersInterface& getGaugeObservablesParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& getWilsonTwoFlavourChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::StaggeredChiralCondensateParametersInterface& getStaggeredChiralCondensateParametersInterface() = 0;
            virtual const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsParametersInterface() = 0;
            virtual const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& getStaggeredTwoFlavourCorrelatorsParametersInterface() = 0;

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

            virtual const physics::SourcesParametersInterface& getSourcesParametersInterface() = 0;

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
            virtual const physics::AdditionalParameters& getWilsonAdditionalParameters(bool) = 0;
            virtual const physics::AdditionalParameters& getStaggeredAdditionalParameters() = 0;
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

    template<> inline const physics::AdditionalParameters& InterfacesHandler::getAdditionalParameters<physics::lattices::Spinorfield>(bool withMassPreconditioning)
    {
        return getWilsonAdditionalParameters(withMassPreconditioning);
    }
    template<> inline const physics::AdditionalParameters& InterfacesHandler::getAdditionalParameters<physics::lattices::Spinorfield_eo>(bool withMassPreconditioning)
    {
        return getWilsonAdditionalParameters(withMassPreconditioning);
    }
    template<> inline const physics::AdditionalParameters& InterfacesHandler::getAdditionalParameters<physics::lattices::Staggeredfield_eo>(bool)
    {
        return getStaggeredAdditionalParameters();
    }
    template<> inline const physics::AdditionalParameters& InterfacesHandler::getAdditionalParameters<physics::lattices::Rooted_Staggeredfield_eo>(bool)
    {
        return getStaggeredAdditionalParameters();
    }

}


