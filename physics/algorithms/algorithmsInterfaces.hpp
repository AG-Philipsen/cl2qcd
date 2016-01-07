/** @file
 * algorithms interfaces declaration
 *
 * Copyright 2016 Alessandro Sciarra, Christopher Czaban
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

#include "../../common_header_files/types.h"

namespace physics{
	namespace algorithms{

	    class SolversParametersInterface {
	        public:
	            virtual ~SolversParametersInterface() {}
	            virtual unsigned getCgMax() const = 0;
	            virtual unsigned getIterRefresh() const = 0;
	            virtual common::solver getSolver() const = 0;
	            virtual unsigned getCgIterationBlockSize() const = 0;
	            virtual unsigned getCgMinimumIterationCount() const = 0;
	            virtual bool getCgUseAsyncCopy() const = 0;
	            virtual bool getUseMergeKernelsSpinor() const = 0;
	    };

        class ForcesParametersInterface{
            public:
                virtual ~ForcesParametersInterface(){}
                virtual common::action getFermact() const = 0;
                virtual double getForcePreconditioning() const = 0;
                virtual double getKappa() const = 0;
                virtual double getKappaMp() const = 0;
                virtual hmc_float getMubar() const = 0;
                virtual hmc_float getMubarMp() const = 0;
                virtual unsigned getRhoIterations() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual bool getUseSmearing() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseRectangles() const = 0;
        };

        class MinMaxEigenvalueParametersInterface {
            public:
                virtual ~MinMaxEigenvalueParametersInterface(){}
                virtual unsigned getFindMinMaxIterationBlockSize() const = 0;
                virtual unsigned getFindMinMaxMaxValue() const = 0;
        };

		class InversionParemetersInterface{
            public:
                virtual ~InversionParemetersInterface(){}
                virtual common::action getFermact() const = 0;
                virtual double getKappa() const = 0;
                virtual double getMubar() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual double getSolverPrec() const = 0;
                virtual bool getUseEo() const = 0;
                virtual bool getUseSmearing() const = 0;
		};

		class IntegratorParametersInterface{
            public:
                virtual ~IntegratorParametersInterface(){}
                virtual unsigned getIntegrationSteps(const unsigned) const = 0;
                virtual common::integrator getIntegrator(const unsigned) const = 0;
                virtual double getKappaMp() const = 0;
                virtual double getLambda(const unsigned) const = 0;
                virtual double getMubarMp() const = 0;
                virtual unsigned getNumTimescales() const = 0;
                virtual double getTau() const = 0;
                virtual bool getUseMp() const = 0;
		};

        class MolecularDynamicsInterface{
            public:
                virtual ~MolecularDynamicsInterface(){}
                virtual double getKappaMp() const = 0;
                virtual double getMubarMp() const = 0;
                virtual double getSolverPrec() const = 0;
        };

		class MetropolisParametersInterface{
            public:
                virtual ~MetropolisParametersInterface(){}
                virtual double getC0() const = 0;
                virtual double getC1() const = 0;
                virtual common::action getFermact() const = 0;
                virtual double getKappa() const = 0;
                virtual double getKappaMp() const = 0;
                virtual double getMubar() const = 0;
                virtual double getMubarMp() const = 0;
                virtual size_t getRectanglesNormalization() const = 0;
                virtual common::solver getSolver() const = 0;
                virtual double getSolverPrec() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
                virtual bool getUseRectangles() const = 0;
		};

		class HmcParametersInterface{
            public:
                virtual ~HmcParametersInterface(){}
                virtual double getBeta() const = 0;
                virtual double getKappa() const = 0;
                virtual double getKappaMp() const = 0;
                virtual double getMubar() const = 0;
                virtual double getMubarMp() const = 0;
                virtual bool getUseEo() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
		};

		class RhmcParametersInterface{
            public:
                virtual ~RhmcParametersInterface(){}
                virtual double getBeta() const = 0;
                virtual bool getConservative() const = 0;
                virtual double getFindMinMaxPrec() const = 0;
                virtual double getKappaMp()	const = 0;
                virtual double getMass() const = 0;
                virtual double getMubarMp() const = 0;
                virtual bool getUseGaugeOnly() const = 0;
                virtual bool getUseMp() const = 0;
                virtual bool getUseEo() const = 0;
		};

	}
}



#include "../../meta/inputparameters.hpp"
#include "../../meta/util.hpp"

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
                virtual double getKappa() const override
                {
                    return parameters.get_kappa();
                }
                virtual double getKappaMp() const override
                {
                    return parameters.get_kappa_mp();
                }
                virtual hmc_float getMubar() const override
                {
                    return meta::get_mubar(parameters);
                }
                virtual hmc_float getMubarMp() const override
                {
                    return meta::get_mubar_mp(parameters);
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
				virtual double getKappa() const override
				{
					return parameters.get_kappa();
				}
				virtual double getMubar() const override
				{
					return meta::get_mubar(parameters);
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
                virtual double getKappaMp() const override
				{
					return parameters.get_kappa_mp();
				}
                virtual double getLambda(const unsigned timescale) const override
                {
                    return parameters.get_lambda(timescale);
                }
				virtual double getMubarMp() const override
				{
					return meta::get_mubar_mp(parameters);
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
            virtual double getKappaMp() const override
            {
                return parameters.get_kappa_mp();
            }
            virtual double getMubarMp() const override
            {
                return meta::get_mubar_mp(parameters);
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
			virtual double getKappa() const override
			{
				return parameters.get_kappa();
			}
			virtual double getKappaMp() const override
			{
				return parameters.get_kappa_mp();
			}
			virtual double getMubar() const override
			{
				return meta::get_mubar(parameters);
			}
			virtual double getMubarMp() const override
			{
				return meta::get_mubar_mp(parameters);
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
                virtual double getKappa() const override
                {
                    return parameters.get_kappa();
                }
                virtual double getKappaMp() const override
                {
                    return parameters.get_kappa_mp();
                }
                virtual double getMubar() const override
                {
                    return meta::get_mubar(parameters);
                }
                virtual double getMubarMp() const override
                {
                    return meta::get_mubar_mp(parameters);
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
			virtual double getKappaMp() const override
			{
				return parameters.get_kappa_mp();
			}
			virtual double getMass() const override
			{
				return parameters.get_mass();
			}
			virtual double getMubarMp() const override
			{
				return meta::get_mubar_mp(parameters);
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
