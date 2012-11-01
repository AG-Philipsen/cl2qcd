/** @file
 * Implementation of the physics::algorithms::solver class
 *
 * (c) 2012 Christopher Pinke <pinke@th.uni-frankfurt.de>
 */

#include "solver.hpp"

#include "../../hardware/system.hpp"
#include "../../logger.hpp"

#include "solver.hpp"

int physics::algorithms::solver::Solver::get_iter_max() const noexcept
{
  return iter_max;
}

int physics::algorithms::solver::Solver::get_iter_refresh() const noexcept
{
  return iter_refresh;
}

int physics::algorithms::solver::Solver::get_iter() 
{
  return iter;
}

hmc_float physics::algorithms::solver::Solver::get_prec() const noexcept
{
  return acc_prec;
}
