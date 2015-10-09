/** @file
 * Interfaces for physics::observables
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

#include "../../meta/inputparameters.hpp"

namespace physics {
    namespace observables {

        class GaugeObservablesParametersInterface {
            public:
                virtual ~GaugeObservablesParametersInterface() {};
                virtual bool measureRectangles() const = 0;
                virtual bool measureTransportCoefficientKappa() const = 0;
                virtual bool printToScreen() const = 0;
                virtual hmc_float getBeta() const = 0;
                virtual std::string getTransportCoefficientKappaFilename() const = 0;
                virtual std::string getRectanglesFilename() const = 0;
                virtual std::string getGaugeObservablesFilename(std::string) const = 0;
                virtual hmc_float getTemporalPlaquetteNormalization() const = 0;
                virtual hmc_float getSpatialPlaquetteNormalization() const = 0;
                virtual hmc_float getPlaquetteNormalization() const = 0;
                virtual unsigned getSpatialVolume() const = 0;
                virtual hmc_float getPolyakovLoopNormalization() const = 0;
        };


        class GaugeObservablesParametersImplementation final : public GaugeObservablesParametersInterface {
            public:
                GaugeObservablesParametersImplementation() = delete;
                GaugeObservablesParametersImplementation(const meta::Inputparameters& parametersIn) : parameters(parametersIn) {}
                ~GaugeObservablesParametersImplementation() {}
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
                hmc_float getTemporalPlaquetteNormalization() const override
                {
                    return meta::get_tplaq_norm(parameters);
                }
                hmc_float getSpatialPlaquetteNormalization() const override
                {
                    return meta::get_splaq_norm(parameters);
                }
                hmc_float getPlaquetteNormalization() const override
                {
                    return meta::get_plaq_norm(parameters);
                }
                unsigned getSpatialVolume() const override
                {
                    return meta::get_volspace(parameters.get_nspace());
                }
                hmc_float getPolyakovLoopNormalization() const override
                {
                    return meta::get_poly_norm(parameters);
                }
            private:
                const meta::Inputparameters& parameters;
        };

    }
}
