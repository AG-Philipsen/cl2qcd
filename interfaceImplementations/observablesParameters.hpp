/** @file
 * observablesParameters.hpp
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

#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"
#include "../physics/observables/observablesInterfaces.hpp"

namespace physics {
    namespace observables {

        class GaugeObservablesParametersImplementation final : public GaugeObservablesParametersInterface {
            public:
                GaugeObservablesParametersImplementation() = delete;
                GaugeObservablesParametersImplementation(const meta::Inputparameters& parametersIn)
                        : parameters(parametersIn)
                {
                }
                ~GaugeObservablesParametersImplementation()
                {
                }
                bool measureRectangles() const
                {
                    return parameters.get_measure_rectangles();
                }
                bool measureTransportCoefficientKappa() const override
                {
                    return parameters.get_measure_transportcoefficient_kappa();
                }
                bool printToScreen() const override
                {
                    return parameters.get_print_to_screen();
                }
                hmc_float getBeta() const override
                {
                    return parameters.get_beta();
                }
                std::string getTransportCoefficientKappaFilename() const override
                {
                    return parameters.get_transportcoefficientKappaFilename();
                }
                std::string getRectanglesFilename() const override
                {
                    return parameters.get_rectanglesFilename();
                }
                std::string getGaugeObservablesFilename(std::string configurationName) const override
                {
                    return meta::get_gauge_obs_file_name(parameters, configurationName);
                }
                unsigned getTemporalPlaquetteNormalization() const override
                {
                    return meta::get_tplaq_norm(parameters);
                }
                unsigned getSpatialPlaquetteNormalization() const override
                {
                    return meta::get_splaq_norm(parameters);
                }
                unsigned getPlaquetteNormalization() const override
                {
                    return meta::get_plaq_norm(parameters);
                }
                unsigned getSpatialVolume() const override
                {
                    return meta::get_volspace(parameters.get_nspace());
                }
                unsigned getPolyakovLoopNormalization() const override
                {
                    return meta::get_poly_norm(parameters);
                }
            private:
                const meta::Inputparameters& parameters;
        };

        class WilsonTwoFlavourChiralCondensateParametersImplementation final : public WilsonTwoFlavourChiralCondensateParametersInterface {
            public:
                WilsonTwoFlavourChiralCondensateParametersImplementation() = delete;
                WilsonTwoFlavourChiralCondensateParametersImplementation(const meta::Inputparameters& parametersIn)
                        : parameters(parametersIn)
                {
                }
                ~WilsonTwoFlavourChiralCondensateParametersImplementation()
                {
                }
                common::action getFermionicActionType() const override
                {
                    return parameters.get_fermact();
                }
                common::pbp_version getPbpVersion() const override
                {
                    return parameters.get_pbp_version();
                }
                bool measurePbp() const override
                {
                    return parameters.get_measure_pbp();
                }
                int getNumberOfSources() const override
                {
                    return parameters.get_num_sources();
                }
                hmc_float getKappa() const override
                {
                    return parameters.get_kappa();
                }
                hmc_float getMubar() const override
                {
                    return meta::get_mubar(parameters);
                }
                unsigned get4dVolume() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace());
                }
                bool useEvenOdd() const override
                {
                    return parameters.get_use_eo();
                }
                std::string getPbpFilename(std::string configurationName) const override
                {
                    return meta::get_ferm_obs_pbp_file_name(parameters, configurationName);
                }

            private:
                const meta::Inputparameters& parameters;

        };

        class StaggeredChiralCondensateParametersImplementation final : public StaggeredChiralCondensateParametersInterface {
            public:
                StaggeredChiralCondensateParametersImplementation() = delete;
                StaggeredChiralCondensateParametersImplementation(const meta::Inputparameters& parametersIn)
                        : parameters(parametersIn)
                {
                }
                ~StaggeredChiralCondensateParametersImplementation()
                {
                }
                bool measurePbp() const override
                {
                    return parameters.get_measure_pbp();
                }
                int getNumberOfSources() const override
                {
                    return parameters.get_num_sources();
                }
                hmc_float getMass() const override
                {
                    return parameters.get_mass();
                }
                hmc_float getSolverPrecision() const override
                {
                    return parameters.get_solver_prec();
                }
                hmc_float getNumberOfTastes() const override
                {
                    return parameters.get_num_tastes();
                }
                unsigned get4dVolume() const override
                {
                    return meta::get_vol4d(parameters.get_ntime(), parameters.get_nspace());
                }
                std::string getPbpFilename(std::string configurationName) const override
                {
                    return meta::get_ferm_obs_pbp_file_name(parameters, configurationName);
                }
                unsigned getPbpNumberOfMeasurements() const override
                {
                    return parameters.get_pbp_measurements();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class WilsonTwoFlavourCorrelatorsParametersImplementation final : public WilsonTwoFlavourCorrelatorsParametersInterface {
            public:
                WilsonTwoFlavourCorrelatorsParametersImplementation() = delete;
                WilsonTwoFlavourCorrelatorsParametersImplementation(const meta::Inputparameters& parametersIn)
                        : parameters(parametersIn)
                {
                }
                ~WilsonTwoFlavourCorrelatorsParametersImplementation()
                {
                }
                bool printToScreen() const override
                {
                    return parameters.get_print_to_screen();
                }
                void printInformationOfFlavourDoubletCorrelator(std::ostream* of = nullptr) const override
                {
                    if(of == nullptr)
                        meta::print_info_flavour_doublet_correlators(parameters);
                    else
                        meta::print_info_flavour_doublet_correlators(of, parameters);
                }
                unsigned getCorrelatorDirection() const override
                {
                    return parameters.get_corr_dir();
                }
                common::sourcetypes getSourceType() const override
                {
                    return parameters.get_sourcetype();
                }
                unsigned getNs() const override
                {
                    return parameters.get_nspace();
                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                std::string getCorrelatorFilename(std::string currentConfigurationName) const override
                {
                    return meta::get_ferm_obs_corr_file_name(parameters, currentConfigurationName);
                }
                bool placeSourcesOnHost() const override
                {
                    return parameters.get_place_sources_on_host();
                }
                int getNumberOfSources() const override
                {
                    return parameters.get_num_sources();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class StaggeredTwoFlavourCorrelatorsParametersImplementation final : public StaggeredTwoFlavourCorrelatorsParametersInterface {
            public:
                StaggeredTwoFlavourCorrelatorsParametersImplementation() = delete;
                StaggeredTwoFlavourCorrelatorsParametersImplementation(const meta::Inputparameters& parametersIn)
                        : parameters(parametersIn)
                {
                }
                ~StaggeredTwoFlavourCorrelatorsParametersImplementation()
                {
                }
                bool printToScreen() const override
                {
                    return parameters.get_print_to_screen();
                }
                void printInformationOfFlavourDoubletCorrelator(std::ostream* of = nullptr) const override
                {
                    if(of == nullptr)
                        meta::print_info_flavour_doublet_correlators(parameters);
                    else
                        meta::print_info_flavour_doublet_correlators(of, parameters);
                }
//                unsigned getCorrelatorDirection() const override
//                {
//                    return parameters.get_corr_dir();
//                }
                common::sourcetypes getSourceType() const override
                  {
                     return parameters.get_sourcetype();
                  }
//                unsigned getNs() const override
//                {
//                    return parameters.get_nspace();
//                }
                unsigned getNt() const override
                {
                    return parameters.get_ntime();
                }
                std::string getCorrelatorFilename(std::string currentConfigurationName) const override
                {
                    return meta::get_ferm_obs_corr_file_name(parameters, currentConfigurationName);
                }
                hmc_float getSolverPrecision() const override
                {
                    return parameters.get_solver_prec();
                }
//                bool placeSourcesOnHost() const override
//                {
//                    return parameters.get_place_sources_on_host();
//                }
                int getNumberOfSources() const override
                {
                    return parameters.get_num_sources();
                }

            private:
                const meta::Inputparameters& parameters;
        };

    }
}

