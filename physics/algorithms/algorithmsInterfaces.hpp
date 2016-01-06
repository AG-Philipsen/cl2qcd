/*
 * algorithmsInterfaces.hpp
 *
 *  Created on: 24 Nov 2015
 *      Author: czaban
 */

#pragma once

#include "../../common_header_files/types.h"

namespace physics{
	namespace algorithms{

		class InversionParemetersInterface{
		public:
			~InversionParemetersInterface(){}
			virtual common::action getFermact() const = 0;
			virtual double getKappa() const = 0;
			virtual double getMubar() const = 0;
			virtual common::solver getSolver() = 0;
			virtual double getSolverPrec() = 0;
			virtual bool getUseEo() = 0;
			virtual bool getUseSmearing() = 0;
		};

		class IntegratorParametersInterface{
			~IntegratorParametersInterface(){}
			virtual unsigned getIntegrationSteps() const = 0;
			virtual common::integrator getIntegrator() const = 0;
			virtual double getKappaMp() const = 0;
			virtual double getLambda() const = 0;
			virtual double getMubarMp() const = 0;
			virtual unsigned getNumTimescales() const = 0;
			virtual double getTau() const = 0;
			virtual bool getUseMp() const = 0;
		};

		class MetropolisParametersInterface{
		public:
			~MetropolisParametersInterface(){}
			virtual double getC0() const = 0;
			virtual double getC1() const = 0;
			virtual common::action getFermact() const = 0;
			virtual double getKappa() const = 0;
			virtual double getKappaMp() const = 0;
			virtual double getMubar() const = 0;
			virtual double getMubarMp() const = 0;
			virtual size_t getRectNorm() const = 0;
			virtual common::solver getSolver() const = 0;
			virtual double getSolverPrec() const = 0;
			virtual bool getUseGaugeOnly() const = 0;
			virtual bool getUseMp() const = 0;
			virtual bool getUseRectangles() const = 0;
		};

		class HmcParametersInterface{
		public:
			~HmcParametersInterface(){}
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
			~RhmcParametersInterface(){}
			virtual double getBeta() const = 0;
			virtual bool getConservative() const = 0;
			virtual double getFindMinMaxPrec() const = 0;
			virtual double getKappaMp()	const = 0;
			virtual double getMass() const = 0;
			virtual double getMubarMp() const = 0;
			virtual bool getUseGaugeOnly() const = 0;
			virtual bool getUseMp() const = 0;
		};

		class MolecularDynamicsInterface{
		public:
			~MolecularDynamicsInterface(){}
			virtual double getKappaMp() const = 0;
			virtual double getMubarMp() const = 0;
			virtual double getSolverPrec() const = 0;
		};
	}
}


namespace physics{
	namespace algorithms{

		class InversionParametersImplementation: public InversionParemetersInterface{
			public:
				~InversionParametersImplementation()
				{
				}
				InversionParametersImplementation() = delete;
				InversionParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
				{
				}
				virtual common::action get_fermact() const override
				{
					return parameters.get_fermact();
				}
				virtual double getKappa() const override
				{
					return paramters.get_kappa();
				}
				virtual double getMubar() const override
				{
					return parameters.get_mubar();
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

		class IntegratorParametersImplementation: public IntegratorParametersInterface{
			public:
				~IntegratorParametersImplementation()
				{
				}
				IntegratorParametersImplementation() = delete;
				IntegratorParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
				{
				}
				virtual unsigned getIntegrationSteps() const = 0;
				virtual common::integrator getIntegrator() const override
				{
					return parameters.get_integrationsteps();
				}
				virtual double getKappaMp() const override
				{
					return parameters.get_integrator();
				}
				virtual double getLambda() const override
				{
					return parameters.get_kappa_mp();
				}
				virtual double getMubarMp() const override
				{
					return meta::get_mubar_mp(parameters)
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
				const meta::Inputparameters& parameters
		};

		class MetropolisParametersImplementation: public MetropolisParametersInterface{
		public:
			~MetropolisParametersImplementation(){}
			MetropolisParametersImplementation() = delete;
			MetropolisParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
			{
			}
			virtual double getC0()
			{
				return meta::get_c0(parameters);
			}
			virtual double getC1()
			{
				return meta::get_c1(parameters);
			}
			virtual common::action getFermact()
			{
				return parameters.get_fermact();
			}
			virtual double getKappa()
			{
				return parameters.get_kappa();
			}
			virtual double getKappaMp()
			{
				return parameters.get_kappa_mp();
			}
			virtual double getMubar()
			{
				return meta::get_mubar(parameters);
			}
			virtual double getMubarMp()
			{
				return meta::get_mubar_mp(parameters);
			}
			virtual size_t getRectNorm()
			{
				return meta::get_rect_norm(parameters);
			}
			virtual common::solver getSolver()
			{
				return parameters.get_solver();
			}
			virtual double getSolverPrec()
			{
				return parameters.get_solver_prec();
			}
			virtual bool getUseGaugeOnly()
			{
				return parameters.get_use_gauge_only();
			}
			virtual bool getUseMp()
			{
				return parameters.get_use_mp();
			}
			virtual bool getUseRectangles()
			{
				return meta::get_use_rectangles();
			}


		private:
			const meta::Inputparameters& parameters
		};

		class HmcParametersImplementation: public HmcParametersInterface{
				public:
					~HmcParametersImplementation(){}
					HmcParametersImplementation() = delete;
					HmcParametersImplementation(meta::Inputparameters& paramsIn): parameters(paramsIn)
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


		class RhmcParametersImplementation: public RhmcParametersInterface{
		public:
			~RhmcParametersImplementation(){}
			RhmcParametersImplementation() = delete;
			RhmcParametersImplementation(const meta::Inputparameters& paramsIn): paramsIn(parameters)
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
				return meta::get_mubar_mp();
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


		class MolecularDynamicsImplementation: public MolecularDynamicsInterface{
		public:
			~MolecularDynamicsImplementation(){}
			MolecularDynamicsImplementation() = delete;
			MolecularDynamicsImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn){}
			virtual double getKappaMp() const override
			{
				return parameters.get_kappa_mp();
			}
			virtual double getMubarMp() const override
			{
				return meta::get_mubar_mp();
			}
			virtual double getSolverPrec() const override
			{
				return parameters.get_solver_prec();
			}

		private:
			const meta::Inputparameters& parameters;
		};
	}
}
