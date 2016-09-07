/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#ifndef _META_PARAMETERS_SOLVER_HPP_
#define _META_PARAMETERS_SOLVER_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersSolver {
public:

	common::solver get_solver() const noexcept;
	common::solver get_solver_mp() const noexcept;
	double get_solver_prec() const noexcept;
	double get_force_prec() const noexcept;
	int get_iter_refresh() const noexcept;
	int get_iter_refresh_mp() const noexcept;
	int get_cg_iteration_block_size() const noexcept;
	bool get_cg_use_async_copy() const noexcept;
	int get_cg_minimum_iteration_count() const noexcept;
	int get_cgmax() const noexcept;
	int get_cgmax_mp() const noexcept;
	bool get_profile_solver() const noexcept;

private:
	po::options_description options;

	double solver_prec;
	double force_prec;
	int iter_refresh;
	int iter_refresh_mp;
	int cgmax;
	int cgmax_mp;
	int cg_iteration_block_size;
	bool cg_use_async_copy;
	int cg_minimum_iteration_count;
	bool profile_solver;

protected:
	ParametersSolver();
	virtual ~ParametersSolver();
	ParametersSolver(ParametersSolver const&) = delete;
	ParametersSolver & operator=(ParametersSolver const&) = delete;
	po::options_description & getOptions();

	//at the moment, only 2 solvers are implemented..
	common::solver _solver;
	common::solver _solver_mp;
};

}

#endif
