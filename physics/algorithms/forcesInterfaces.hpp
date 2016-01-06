/*
 * ForcesParametersInterface.hpp
 *
 *  Created on: 23 Nov 2015
 *      Author: czaban
 */

#pragma once

#include "../../common_header_files/types.h"

namespace physics{
	namespace algorithms {

		class ForcesParametersInterfaces{
			public:
				virtual ~ForcesParametersInterfaces(){}
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


		class ForcesParametersImplementation: public ForcesParametersInterfaces{
			public:
				ForcesParametersImplementation() = delete;
				ForcesParametersImplementation(const meta::Inputparameters& paramsIn): parameters(paramsIn)
				{
				}
				~ForcesParametersImplementation()
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
					return meta::get_mubar_mp(parameters)
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
					return paramters.get_use_rectangles();
				}

			private:
				const meta::Inputparameters& parameters;
		};


	}
}
