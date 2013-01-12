/** @file
 * Declaration of the physics::lattices::Spinorfield_eo class
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_EO_
#define _PHYSICS_LATTICES_SPINORFIELD_EO_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/spinor.hpp"
#include "../prng.hpp"
#include "spinorfield.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a gaugefield.
 */
class Spinorfield_eo {

public:
	/**
	 * Construct a gaugefield based on the input-files of the system
	 */
	Spinorfield_eo(const hardware::System&);

	/**
	 * Release resources
	 */
	~Spinorfield_eo();

	/*
	 * Spinorfield_eos cannot be copied
	 */
	Spinorfield_eo& operator=(const Spinorfield_eo&) = delete;
	Spinorfield_eo(const Spinorfield_eo&) = delete;
	Spinorfield_eo() = delete;

	/**
	 * Get the buffers containing the gaugefield state on the devices.
	 */
	const std::vector<const hardware::buffers::Spinor *> get_buffers() const noexcept;

	/**
	 * Apply Gamma5 to the Spinorfield
	 */
	void gamma5() const;

	/**
	 * Set Spinorfield to zero
	 */
	void zero() const;

	/**
	 * Set Spinorfield to cold
	 */
	void cold() const;

	/**
	 * Set Spinorfield to be gaussian.
	 */
	void gaussian(const physics::PRNG& prng) const;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Spinor *> buffers;
};

/**
 * Calculate the scalar product of two spinorfields.
 */
hmc_complex scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right);

/**
 * Calculate the squarenorm of the spinorfield
 */
hmc_float squarenorm(const Spinorfield_eo& field);

/**
 * Perform the BLAS operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const Spinorfield_eo& y);

/**
 * Perform the BLAS operation sax.
 *
 * out = alpha * x
 */
void sax(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x);

/**
 * Perform the BLAS operation saxsbypz.
 *
 * out = alpha * x + beta * y + z
 */
void saxsbypz(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const hmc_complex beta, const Spinorfield_eo& y, const Spinorfield_eo& z);

/**
 * Split the given Spinorfield into even and odd Spinorfield_eo.
 *
 * @param[out] even The even part
 * @param[out] odd  The odd part
 * @param[in]  in   The Spinorfield to split
 */
void convert_to_eoprec(const Spinorfield_eo* even, const Spinorfield_eo* odd, const Spinorfield& in);

/**
 * Merge the given even and odd Spinorfield_eo into one Spinorfield.
 *
 * @param[out] merged The merged Spinorfield
 * @param[in]  even   The even part
 * @param[in]  odd    The odd part
 */
void convert_from_eoprec(const Spinorfield* merged, const Spinorfield_eo& even, const Spinorfield_eo& odd);
}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_EO_ */
