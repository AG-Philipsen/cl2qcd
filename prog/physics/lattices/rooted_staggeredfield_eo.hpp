/** @file
 * Declaration of the physics::lattices::Rooted_Staggeredfield_eo class
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#ifndef _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_
#define _PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_

#include "../../hardware/system.hpp"
#include "../algorithms/rational_approximation.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a rooted staggeredfield (with eo preconditioning).
 */
class Rooted_Staggeredfield_eo : public Staggeredfield_eo, public physics::algorithms::Rational_Coefficients{

public:
	/**
	 * Construct a rooted staggeredfield based on the input-files of the system
	 */
	Rooted_Staggeredfield_eo(const hardware::System&);

	/**
	 * Staggeredfield_eo cannot be copied
	 */
	Rooted_Staggeredfield_eo& operator=(const Rooted_Staggeredfield_eo&) = delete;
	Rooted_Staggeredfield_eo(const Rooted_Staggeredfield_eo&) = delete;
	Rooted_Staggeredfield_eo() = delete;

	/**
	 * Rescale coefficients on the basis of a Rational_Approximation objects
	 */
	void Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const hardware::System& system, hmc_float prec, bool conservative=false);

};

}
}

#endif /*_PHYSICS_LATTICES_ROOTED_STAGGEREDFIELD_EO_ */

