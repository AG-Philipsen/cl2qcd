/** @file
 * interfacesHandler.hpp
 *
 * Copyright 2016 Alessandro Sciarra
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

#include <memory>
#include "physicsParameters.hpp"
#include "algorithmsParameters.hpp"
#include "latticesParameters.hpp"
#include "fermionmatrixParameters.hpp"
#include "observablesParameters.hpp"
#include "../physics/interfacesHandler.hpp"

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
                  staggeredTwoFlavourCorrelatorsParametersInterface{nullptr},
                  solversParametersInterface{nullptr},
                  minMaxEigenvalueParametersInterface{nullptr},
                  forcesParametersInterface{nullptr},
                  inversionParametersInterface{nullptr},
                  integratorParametersInterface{nullptr},
                  molecularDynamicsParametersInterface{nullptr},
                  metropolisParametersInterface{nullptr},
                  hmcParametersInterface{nullptr},
                  rhmcParametersInterface{nullptr},
                  sourcesParametersInterface{nullptr},
                  wilsonAdditionalParameters{nullptr},
                  wilsonAdditionalParametersMp{nullptr},
                  staggeredAdditionalParameters{nullptr} {}
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
            const physics::observables::WilsonTwoFlavourCorrelatorsParametersInterface& getWilsonTwoFlavourCorrelatorsParametersInterface() override
            {
                if(wilsonTwoFlavourCorrelatorsParametersInterface == nullptr)
                    wilsonTwoFlavourCorrelatorsParametersInterface = std::unique_ptr<const physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation>(new physics::observables::WilsonTwoFlavourCorrelatorsParametersImplementation{parameters});
                return *wilsonTwoFlavourCorrelatorsParametersInterface;
            }
            const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface& getStaggeredTwoFlavourCorrelatorsParametersInterface() override
            {
                if(staggeredTwoFlavourCorrelatorsParametersInterface == nullptr)
                    staggeredTwoFlavourCorrelatorsParametersInterface = std::unique_ptr<const physics::observables::StaggeredTwoFlavourCorrelatorsParametersImplementation>(new physics::observables::StaggeredTwoFlavourCorrelatorsParametersImplementation{parameters});
                return *staggeredTwoFlavourCorrelatorsParametersInterface;
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

            const physics::SourcesParametersInterface& getSourcesParametersInterface()
            {
                if( sourcesParametersInterface == nullptr)
                    sourcesParametersInterface = std::unique_ptr<const physics::SourcesParametersImplementation>(new physics::SourcesParametersImplementation{parameters});
                return *sourcesParametersInterface;
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

            virtual const physics::AdditionalParameters& getWilsonAdditionalParameters(bool withMassPreconditioning) override
            {
                if(withMassPreconditioning){
                    if(wilsonAdditionalParametersMp == nullptr)
                        wilsonAdditionalParametersMp = std::unique_ptr<const physics::WilsonAdditionalParameters>(new physics::WilsonAdditionalParameters{parameters, true});
                    return *wilsonAdditionalParametersMp;
                }else{
                    if(wilsonAdditionalParameters == nullptr)
                        wilsonAdditionalParameters = std::unique_ptr<const physics::WilsonAdditionalParameters>(new physics::WilsonAdditionalParameters{parameters, false});
                    return *wilsonAdditionalParameters;
                }
            }

            virtual const physics::AdditionalParameters& getStaggeredAdditionalParameters() override
            {
                if(staggeredAdditionalParameters == nullptr)
                    staggeredAdditionalParameters = std::unique_ptr<const physics::StaggeredAdditionalParameters>(new physics::StaggeredAdditionalParameters{parameters});
                return *staggeredAdditionalParameters;
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
            std::unique_ptr<const physics::observables::StaggeredTwoFlavourCorrelatorsParametersInterface> staggeredTwoFlavourCorrelatorsParametersInterface;
            std::unique_ptr<const physics::algorithms::SolversParametersInterface> solversParametersInterface;
            std::unique_ptr<const physics::algorithms::MinMaxEigenvalueParametersInterface> minMaxEigenvalueParametersInterface;
            std::unique_ptr<const physics::algorithms::ForcesParametersInterface> forcesParametersInterface;
            std::unique_ptr<const physics::algorithms::InversionParemetersInterface> inversionParametersInterface;
            std::unique_ptr<const physics::algorithms::IntegratorParametersInterface> integratorParametersInterface;
            std::unique_ptr<const physics::algorithms::MolecularDynamicsInterface> molecularDynamicsParametersInterface;
            std::unique_ptr<const physics::algorithms::MetropolisParametersInterface> metropolisParametersInterface;
            std::unique_ptr<const physics::algorithms::HmcParametersInterface> hmcParametersInterface;
            std::unique_ptr<const physics::algorithms::RhmcParametersInterface> rhmcParametersInterface;
            std::unique_ptr<const physics::SourcesParametersInterface> sourcesParametersInterface;
            std::unique_ptr<const physics::AdditionalParameters> wilsonAdditionalParameters;
            std::unique_ptr<const physics::AdditionalParameters> wilsonAdditionalParametersMp;
            std::unique_ptr<const physics::AdditionalParameters> staggeredAdditionalParameters;

    };

    /*
     * NOTE: In the InterfacesHandlerImplementation we need two wilsonAdditionalParameters pointers, one without mass preconditioning and one with.
     *       The reason for this is that each interface of the InterfacesHandler is built only the first time the corresponding getter is called.
     *       To understand more in detail which is the problem, let us suppose to have only one WilsonAdditionalParameter pointer. This is allocated
     *       when the getter is called for the first time and later, if the getter is called again, the pointed object does not change, but it is just
     *       returned. But the WilsonAdditionalParameter getter has a bool argument to ask either for the mass preconditioned additional parameters
     *       or to ask for the normal additional parameters. So if this boolean variable changes in two successive call to the getter, having only one
     *       pointer would mean to basically ignore this bool variable.
     *       So we must have an object to be returned when this bool is true and one when it is false.
     */
}
