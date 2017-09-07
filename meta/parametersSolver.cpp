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

#include "parametersSolver.hpp"

int meta::ParametersSolver::get_cgmax() const noexcept
{
	return cgmax;
}
int meta::ParametersSolver::get_cgmax_mp() const noexcept
{
	return cgmax_mp;
}

double meta::ParametersSolver::get_solver_prec() const noexcept
{
	return solver_prec;
}
double meta::ParametersSolver::get_force_prec() const noexcept
{
	return force_prec;
}
int meta::ParametersSolver::get_iter_refresh() const noexcept
{
	return iter_refresh;
}
int meta::ParametersSolver::get_iter_refresh_mp() const noexcept
{
	return iter_refresh_mp;
}
common::solver meta::ParametersSolver::get_solver() const noexcept
{
	return _solver;
}
common::solver meta::ParametersSolver::get_solver_mp() const noexcept
{
	return _solver_mp;
}

int meta::ParametersSolver::get_cg_iteration_block_size() const noexcept
{
	return cg_iteration_block_size;
}
bool meta::ParametersSolver::get_cg_use_async_copy() const noexcept
{
	return cg_use_async_copy;
}
int meta::ParametersSolver::get_cg_minimum_iteration_count() const noexcept
{
	return cg_minimum_iteration_count;
}

bool meta::ParametersSolver::get_profile_solver() const noexcept
{
	return profile_solver;
}

meta::ParametersSolver::ParametersSolver()
	: options("Solver options")
{
	options.add_options()
	("solver", po::value<std::string>()->default_value("bicgstab"))
	("solver_mp", po::value<std::string>()->default_value("bicgstab"))
	("cgmax", po::value<int>(&cgmax)->default_value(1000))
	("cgmax_mp", po::value<int>(&cgmax_mp)->default_value(1000))
#ifdef _USEDOUBLEPREC_
	("solver_prec", po::value<double>(&solver_prec)->default_value(1e-23))
	("force_prec", po::value<double>(&force_prec)->default_value(1e-12))
#else
	("solver_prec", po::value<double>(&solver_prec)->default_value(1e-16))
	("force_prec", po::value<double>(&force_prec)->default_value(1e-8))
#endif
	("iter_refresh", po::value<int>(&iter_refresh)->default_value(100))
	("iter_refresh_mp", po::value<int>(&iter_refresh_mp)->default_value(100))
	//this is not used. Remove!
	("profile_solver", po::value<bool>(&profile_solver)->default_value(false))
	("cg_iteration_block_size", po::value<int>(&cg_iteration_block_size)->default_value(10), "CG will check the residual only every N iterations")
	("cg_use_async_copy", po::value<bool>(&cg_use_async_copy)->default_value(false), "CG will use residual of iteration N - block_size for termination condition.")
	("cg_minimum_iteration_count", po::value<int>(&cg_minimum_iteration_count)->default_value(0), "CG will perform at least this many itertions. USE ONLY FOR BENCHMARKS!");
}

meta::ParametersSolver::~ParametersSolver() = default;

po::options_description & meta::ParametersSolver::getOptions()
{
	return options;
}
