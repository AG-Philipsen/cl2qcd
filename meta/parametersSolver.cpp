/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#include "parametersSolver.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

static common::solver translateSolverToEnum(std::string);

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
    return translateSolverToEnum(_solver);
}
common::solver meta::ParametersSolver::get_solver_mp() const noexcept
{
    return translateSolverToEnum(_solver_mp);
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

meta::ParametersSolver::ParametersSolver() : options("Solver options")
{
    // clang-format off
	options.add_options()
	("solver", po::value<std::string>(&_solver)->default_value("bicgstab"),"Which type of (restarted) solver to use (one among 'cg', 'bicgstab' and 'bicgstab_save').")
	("solverMP", po::value<std::string>(&_solver_mp)->default_value("bicgstab"),"Which type of solver to use with Mass Preconditioning (one among 'cg', 'bicgstab' and 'bicgstab_save').")
	("cgmax", po::value<int>(&cgmax)->default_value(1000),"The maximum number of iterations in the solver.")
	("cgmaxMP", po::value<int>(&cgmax_mp)->default_value(1000),"The maximum number of iterations in the solver with Mass Preconditioning.")
#ifdef _USEDOUBLEPREC_
	("solverPrecision", po::value<double>(&solver_prec)->default_value(1e-23),"The precision used in Metropolis inversions.")
	("forcePrecision", po::value<double>(&force_prec)->default_value(1e-12),"The precision used in  Molecular Dynamics inversions.")
#else
	("solverPrecision", po::value<double>(&solver_prec)->default_value(1e-16),"The precision used in Metropolis inversions.")
	("forcePrecision", po::value<double>(&force_prec)->default_value(1e-8),"The precision used in  Molecular Dynamics inversions.")
#endif
	("restartEvery", po::value<int>(&iter_refresh)->default_value(100),"The frequency at which the current approximate solution becomes the new initial guess for the next 'restartEvery' iterations of the solver.")
	("restartEveryMP", po::value<int>(&iter_refresh_mp)->default_value(100),"The frequency at which the current approximate solution becomes the new initial guess for the next 'restartEvery' iterations of the solver, with Mass Preconditioning.")
	("solverResidualCheckEvery", po::value<int>(&cg_iteration_block_size)->default_value(10), "The frequency at which the solver will check the residual.")
	("cgUseAsyncCopy", po::value<bool>(&cg_use_async_copy)->default_value(false), "Whether the solver uses residual of iteration N - 'checkResidualEvery' for termination condition on iteration N.")
	("cgMinimumNumberOfIterations", po::value<int>(&cg_minimum_iteration_count)->default_value(0), "The minimum number of iterations to be performed by the cg solver. To be used for benchmark purposes only!");
    // clang-format on
}

static common::solver translateSolverToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::solver> m;
    m["cg"]            = common::cg;
    m["bicgstab"]      = common::bicgstab;
    m["bicgstab_save"] = common::bicgstab_save;
    common::solver a   = m[s];
    if (a) {
        return a;
    } else {
        throw Invalid_Parameters("Unkown solver!", "cg, bicgstab, bicgstab_save", s);
    }
}
