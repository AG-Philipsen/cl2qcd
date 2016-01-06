/*
 * findMinMaxEigenvalueInterfaces.hpp
 *
 *  Created on: 24 Nov 2015
 *      Author: czaban
 */

#pragma once

#include "../../meta/inputparameters.hpp"

namespace physics {
	namespace algorithms{

		class MinMaxEigenvalueParametersInterface {
			public:
				virtual ~MinMaxEigenvalueInterface{}
				virtual unsigned getFindMinMaxIterationBlockSize() const = 0;
				virtual unsigned getFindMinMaxMaxValue() const = 0;
		};

		class MinMaxEigenvalueParametersImplementation: public MinMaxEigenvalueParametersInterface {
			public:
				virtual ~MinMaxEigenvalueImplementation
				{
				}
				MinMaxEigenvalueImplementation() = delete;
				MinMaxEigenvalueImplementation(const meta::Inputparameters paramsIn): paramters(paramsIn)
				{
				}
				virtual unsigned getFindMinMaxIterationBlockSize() const override
				{
					return parameters.get_findminmax_iteration_block_size();
				}
				virtual unsigned getFindMinMaxValue() const = override
				{
					return paramters.get_findminmax_max();
				}

			private:
				const meta::Inputparameters& parameters;
		};

	}
}
