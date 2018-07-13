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
common::solver meta::ParametersSolver::get_solver() const noexcept
{
    return _solver;
}
common::solver meta::ParametersSolver::get_solver_mp() const noexcept
{
    return _solverMP;
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

meta::ParametersSolver::ParametersSolver()
    :
#ifdef _USEDOUBLEPREC_
    solver_prec(1e-23)
    , force_prec(1e-12)

#else
    solver_prec(1e-16)
    , force_prec(1e-8)
#endif
    , iter_refresh(100)
    , cgmax(1000)
    , cgmax_mp(1000)
    , cg_iteration_block_size(10)
    , cg_use_async_copy(false)
    , cg_minimum_iteration_count(0)
    , options("Solver options")
    , _solverString("bicgstab")
    , _solverMPString("bicgstab")
    , _solver(common::solver::bicgstab)
    , _solverMP(common::solver::bicgstab)
{
    // clang-format off
    options.add_options()
    ("solver", po::value<std::string>(&_solverString)->default_value(_solverString),"Which type of (restarted) solver to use (one among 'cg', 'bicgstab' and 'bicgstab_save').")
    ("solverMP", po::value<std::string>(&_solverMPString)->default_value(_solverMPString),"Which type of solver to use with Mass Preconditioning (one among 'cg', 'bicgstab' and 'bicgstab_save').")
    ("solverMaxIterations", po::value<int>(&cgmax)->default_value(cgmax),"The maximum number of iterations in the solver.")
    ("solverMinIterations", po::value<int>(&cg_minimum_iteration_count)->default_value(cg_minimum_iteration_count), "The minimum number of iterations to be performed by the cg solver. To be used for benchmark purposes only!")
    ("solverMaxIterationsMP", po::value<int>(&cgmax_mp)->default_value(cgmax_mp),"The maximum number of iterations in the solver with Mass Preconditioning.")
    ("solverPrecision", po::value<double>(&solver_prec)->default_value(solver_prec, meta::getDefaultForHelper(solver_prec)),"The precision used in all inversions, except those in the Molecular Dynamics.")
    ("solverForcePrecision", po::value<double>(&force_prec)->default_value(force_prec, meta::getDefaultForHelper(force_prec)),"The precision used in Molecular Dynamics inversions.")
    ("solverRestartEvery", po::value<int>(&iter_refresh)->default_value(iter_refresh),"Every how many iterations the residuum is set to \"A*x-b\" using the current approximate solution before being normally updated.")
    ("solverResiduumCheckEvery", po::value<int>(&cg_iteration_block_size)->default_value(cg_iteration_block_size), "The frequency at which the solver will check the residuum.")
    ("solverUseAsyncCopy", po::value<bool>(&cg_use_async_copy)->default_value(cg_use_async_copy), "Whether the solver uses residuum of iteration N - 'checkResidualEvery' for termination condition on iteration N.");
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

void meta::ParametersSolver::makeNeededTranslations()
{
    _solver   = translateSolverToEnum(_solverString);
    _solverMP = translateSolverToEnum(_solverMPString);
}
