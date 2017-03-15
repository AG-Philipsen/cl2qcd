/** @file
 * Declaration of the physics::lattices::Rooted_Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#ifndef _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_
#define _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_

#include "../../hardware/system.hpp"
#include "../algorithms/rational_approximation.hpp"
//This is to make the template pseudo_randomize friend of this class
#include "util.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a rooted staggeredfield (with eo preconditioning).
 */
class Rooted_Staggeredfield_eo : public Staggeredfield_eo{

public:
	/**
	 * Construct a rooted staggeredfield based on the input-files of the system
	 */
	Rooted_Staggeredfield_eo(const hardware::System&, const RootedStaggeredfieldEoParametersInterface&);
	Rooted_Staggeredfield_eo(const hardware::System&, const RootedStaggeredfieldEoParametersInterface&, const physics::algorithms::Rational_Approximation& approx);

	virtual ~Rooted_Staggeredfield_eo(){}

	/**
	 * Staggeredfield_eo cannot be copied
	 */
	Rooted_Staggeredfield_eo& operator=(const Rooted_Staggeredfield_eo&) = delete;
	Rooted_Staggeredfield_eo(const Rooted_Staggeredfield_eo&) = delete;
	Rooted_Staggeredfield_eo() = delete;

	/**
	 * Rescale coefficients on the basis of a Rational_Approximation objects
	 *  @param A The fermionic operator to calculate the eigenvalues from
	 *  @param gf The gaugefield which A depends on
	 *  @param system The system it is operated on
	 *  @param prec The precision to calculate the eigenvalues up
	 *  @param conservative If true, the maximum eigenvalue found by find_maxmin_eigenvalue is
	 *                      increased by 5% and the minimum one is set to the squared mass of the
	 *                      fermions. This circumvents possible numeric errors.
	 */
	void Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const hmc_float minEigenvalue, const hmc_float maxEigenvalue);

	/**
	 * This method returns the order of the approximation
	 */
	const physics::algorithms::Rational_Coefficients getRationalCoefficients() const noexcept;

private:
	physics::algorithms::Rational_Coefficients rationalCoefficients;
};

}
}

#endif /*_PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_ */

