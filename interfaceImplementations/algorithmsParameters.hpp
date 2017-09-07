/** @file
 * algorithmsParameters.hpp
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
#include "../physics/algorithms/algorithmsInterfaces.hpp"

namespace physics{
    namespace algorithms{

        class SolversParametersImplementation final : public SolversParametersInterface {
            public:
                SolversParametersImplementation() = delete;
                SolversParametersImplementation(const meta::Inputparameters& paramsIn)
                        :parameters(paramsIn)
                {
                }
                virtual ~SolversParametersImplementation()
                {
                }
                virtual unsigned getCgMax() const override
                {
                    return parameters.get_cgmax();
                }
                virtual unsigned getIterRefresh() const override
                {
                    return parameters.get_iter_refresh();
                }
                virtual common::solver getSolver() const override
                {
                    return parameters.get_solver();
                }
                virtual unsigned getCgIterationBlockSize() const override
                {
                    return parameters.get_cg_iteration_block_size();
                }
                virtual unsigned getCgMinimumIterationCount() const override
                {
                    return parameters.get_cg_minimum_iteration_count();
                }
                virtual bool getCgUseAsyncCopy() const override
                {
                    return parameters.get_cg_use_async_copy();
                }
                virtual bool getUseMergeKernelsSpinor() const override
                {
                    return parameters.get_use_merge_kernels_spinor();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class MinMaxEigenvalueParametersImplementation final : public MinMaxEigenvalueParametersInterface {
            public:
                MinMaxEigenvalueParametersImplementation() = delete;
                MinMaxEigenvalueParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
                {
                }
                virtual ~MinMaxEigenvalueParametersImplementation()
                {
                }
                virtual unsigned getFindMinMaxIterationBlockSize() const override
                {
                    return parameters.get_findminmax_iteration_block_size();
                }
                virtual unsigned getFindMinMaxMaxValue() const override
                {
                    return parameters.get_findminmax_max();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class ForcesParametersImplementation final : public ForcesParametersInterface{
            public:
                ForcesParametersImplementation() = delete;
                ForcesParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
                {
                }
                virtual ~ForcesParametersImplementation()
                {
                }
                virtual common::action getFermact() const override
                {
                    return parameters.get_fermact();
                }
                virtual double getForcePreconditioning() const override
                {
                    return parameters.get_force_prec();
                }
                virtual unsigned getRhoIterations() const override
                {
                    return parameters.get_rho_iter();
                }
                virtual common::solver getSolver() const override
                {
                    return parameters.get_solver();
                }
                virtual bool getUseSmearing() const override
                {
                    return parameters.get_use_smearing();
                }
                virtual bool getUseGaugeOnly() const override
                {
                    return parameters.get_use_gauge_only();
                }
                virtual bool getUseRectangles() const override
                {
                    return meta::get_use_rectangles(parameters);
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class InversionParametersImplementation final : public InversionParemetersInterface{
            public:
                InversionParametersImplementation() = delete;
                InversionParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
                {
                }
                ~InversionParametersImplementation()
                {
                }
                virtual common::action getFermact() const override
                {
                    return parameters.get_fermact();
                }
                virtual common::solver getSolver() const override
                {
                    return parameters.get_solver();
                }
                virtual double getSolverPrec() const override
                {
                    return parameters.get_solver_prec();
                }
                virtual bool getUseEo() const override
                {
                    return parameters.get_use_eo();
                }
                virtual bool getUseSmearing() const override
                {
                    return parameters.get_use_smearing();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class IntegratorParametersImplementation final : public IntegratorParametersInterface{
            public:
                IntegratorParametersImplementation() = delete;
                IntegratorParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
                {
                }
                virtual ~IntegratorParametersImplementation()
                {
                }
                virtual unsigned getIntegrationSteps(const unsigned timescale) const override
                {
                    return parameters.get_integrationsteps(timescale);
                }
                virtual common::integrator getIntegrator(const unsigned timescale) const override
                {
                    return parameters.get_integrator(timescale);
                }
                virtual double getLambda(const unsigned timescale) const override
                {
                    return parameters.get_lambda(timescale);
                }
                virtual unsigned getNumTimescales() const override
                {
                    return parameters.get_num_timescales();
                }
                virtual double getTau() const override
                {
                    return parameters.get_tau();
                }
                virtual bool getUseMp() const override
                {
                    return parameters.get_use_mp();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class MolecularDynamicsImplementation final : public MolecularDynamicsInterface{
        public:
            MolecularDynamicsImplementation() = delete;
            MolecularDynamicsImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
            {
            }
            virtual ~MolecularDynamicsImplementation()
            {
            }
            virtual double getSolverPrec() const override
            {
                return parameters.get_solver_prec();
            }

        private:
            const meta::Inputparameters& parameters;
        };

        class MetropolisParametersImplementation final : public MetropolisParametersInterface{
        public:
            MetropolisParametersImplementation() = delete;
            MetropolisParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
            {
            }
            virtual ~MetropolisParametersImplementation()
            {
            }
            virtual double getC0() const override
            {
                return meta::get_c0(parameters);
            }
            virtual double getC1() const override
            {
                return meta::get_c1(parameters);
            }
            virtual common::action getFermact() const override
            {
                return parameters.get_fermact();
            }
            virtual size_t getRectanglesNormalization() const override
            {
                return meta::get_rect_norm(parameters);
            }
            virtual common::solver getSolver() const override
            {
                return parameters.get_solver();
            }
            virtual double getSolverPrec() const override
            {
                return parameters.get_solver_prec();
            }
            virtual bool getUseGaugeOnly() const override
            {
                return parameters.get_use_gauge_only();
            }
            virtual bool getUseMp() const override
            {
                return parameters.get_use_mp();
            }
            virtual bool getUseRectangles() const override
            {
                return meta::get_use_rectangles(parameters);
            }

        private:
            const meta::Inputparameters& parameters;
        };

        class HmcParametersImplementation final : public HmcParametersInterface{
            public:
                HmcParametersImplementation() = delete;
                HmcParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
                {
                }
                virtual ~HmcParametersImplementation()
                {
                }
                virtual double getBeta() const override
                {
                    return parameters.get_beta();
                }
                virtual bool getUseEo() const override
                {
                    return parameters.get_use_eo();
                }
                virtual bool getUseGaugeOnly() const override
                {
                    return parameters.get_use_gauge_only();
                }
                virtual bool getUseMp() const override
                {
                    return parameters.get_use_mp();
                }

            private:
                const meta::Inputparameters& parameters;
        };

        class RhmcParametersImplementation final : public RhmcParametersInterface{
        public:
            RhmcParametersImplementation() = delete;
            RhmcParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
            {
            }
            virtual ~RhmcParametersImplementation()
            {
            }
            virtual double getBeta() const override
            {
                return parameters.get_beta();
            }
            virtual bool getConservative() const override
            {
                return parameters.get_conservative();
            }
            virtual double getFindMinMaxPrec() const override
            {
                return parameters.get_findminmax_prec();
            }
            virtual double getMass() const override
            {
                return parameters.get_mass();
            }
            virtual bool getUseGaugeOnly() const override
            {
                return parameters.get_use_gauge_only();
            }
            virtual bool getUseMp() const override
            {
                return parameters.get_use_mp();
            }
            virtual bool getUseEo() const override
            {
                return parameters.get_use_eo();
            }

        private:
            const meta::Inputparameters& parameters;
        };

    }
}
