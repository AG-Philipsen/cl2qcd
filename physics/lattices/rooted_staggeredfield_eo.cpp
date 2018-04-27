/** @file
 * Implementation of the physics::lattices::Rooted_Staggeredfield_eo class
 *
 * Copyright (c) 2013,2015-2018 Alessandro Sciarra
 * Copyright (c) 2013 Matthias Bach
 * Copyright (c) 2017 Francesca Cuteri
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

#include "rooted_staggeredfield_eo.hpp"

#include <algorithm>  // std::max

physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(
    const hardware::System& system,
    const RootedStaggeredfieldEoParametersInterface& rootedStaggeredfieldEoParametersInterface)
    : rationalCoefficients(
          std::max(rootedStaggeredfieldEoParametersInterface.getMetropolisRationalApproximationOrder(),
                   rootedStaggeredfieldEoParametersInterface.getMolecularDynamicsRationalApproximationOrder()))
{
    const unsigned int numberOfPseudofermions = rootedStaggeredfieldEoParametersInterface.getNumberOfPseudofermions();
    pseudofermions.reserve(numberOfPseudofermions);
    for (unsigned int j = 0; j < numberOfPseudofermions; j++) {
        pseudofermions.emplace_back(
            new physics::lattices::Staggeredfield_eo(system, rootedStaggeredfieldEoParametersInterface));
    }
}

physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(
    const hardware::System& system,
    const RootedStaggeredfieldEoParametersInterface& rootedStaggeredfieldEoParametersInterface,
    const physics::algorithms::Rational_Approximation& approx)
    : rationalCoefficients(approx.Get_order(), approx.Get_a0(), approx.Get_a(), approx.Get_b())
{
    const unsigned int numberOfPseudofermions = rootedStaggeredfieldEoParametersInterface.getNumberOfPseudofermions();
    pseudofermions.reserve(numberOfPseudofermions);
    for (unsigned int j = 0; j < numberOfPseudofermions; j++) {
        pseudofermions.emplace_back(
            new physics::lattices::Staggeredfield_eo(system, rootedStaggeredfieldEoParametersInterface));
    }
}

void physics::lattices::Rooted_Staggeredfield_eo::Rescale_Coefficients(
    const physics::algorithms::Rational_Approximation& approx, const hmc_float minEigenvalue,
    const hmc_float maxEigenvalue)
{
    rationalCoefficients = approx.Rescale_Coefficients(minEigenvalue, maxEigenvalue);
}

unsigned int physics::lattices::Rooted_Staggeredfield_eo::getOrder() const
{
    return rationalCoefficients.Get_order();
}

hmc_float physics::lattices::Rooted_Staggeredfield_eo::get_a0() const
{
    return rationalCoefficients.Get_a0();
}

std::vector<hmc_float> physics::lattices::Rooted_Staggeredfield_eo::get_a() const
{
    return rationalCoefficients.Get_a();
}

std::vector<hmc_float> physics::lattices::Rooted_Staggeredfield_eo::get_b() const
{
    return rationalCoefficients.Get_b();
}

const std::unique_ptr<physics::lattices::Staggeredfield_eo>& physics::lattices::Rooted_Staggeredfield_eo::
operator[](unsigned int index) const
{
    return pseudofermions[index];
}

std::vector<std::unique_ptr<physics::lattices::Staggeredfield_eo>>::iterator
physics::lattices::Rooted_Staggeredfield_eo::begin()
{
    return pseudofermions.begin();
}

std::vector<std::unique_ptr<physics::lattices::Staggeredfield_eo>>::const_iterator
physics::lattices::Rooted_Staggeredfield_eo::begin() const
{
    return pseudofermions.begin();
}

std::vector<std::unique_ptr<physics::lattices::Staggeredfield_eo>>::iterator
physics::lattices::Rooted_Staggeredfield_eo::end()
{
    return pseudofermions.end();
}

std::vector<std::unique_ptr<physics::lattices::Staggeredfield_eo>>::const_iterator
physics::lattices::Rooted_Staggeredfield_eo::end() const
{
    return pseudofermions.end();
}
