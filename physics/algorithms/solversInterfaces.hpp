/*
 * solversInterfaces.hpp
 *
 *  Created on: 23 Nov 2015
 *      Author: czaban
 */

#pragma once

#include "../../common_header_files/types.h"

namespace physics {
	namespace algorithms {

	class SolversParametersInterface {
	public:
		virtual ~SolverParametersInterface() {}
		virtual unsigned getCgMax() const = 0;
		virtual unsigned getIterRefresh() = 0;
		virtual common::solver getSolver() = 0;

		virtual unsigned getCgIterationBlockSize() = 0;
		virtual unsigned getCgMinimumIterationCount() = 0;
		virtual bool getCgUseAsyncCopy() = 0;
		virtual bool getUseMergeKernelsSpinor() = 0;
	};


#include "../../meta/inputparameters.hpp"


	class SolversParametersImplementation: public SolversParametersInterface {
		public:
			SolversParametersImplementation() = delete;
			SolversParametersImplementation(const meta::InputParameters& paramsIn)
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
			const meta::InputParameters& parameters;
	};

	}
}
