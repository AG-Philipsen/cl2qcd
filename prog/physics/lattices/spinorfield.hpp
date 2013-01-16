/** @file
 * Declaration of the physics::lattices::Spinorfield class
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_
#define _PHYSICS_LATTICES_SPINORFIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/plain.hpp"
#include "../prng.hpp"
#include "scalar.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a gaugefield.
 */
class Spinorfield {

public:
	/**
	 * Construct a gaugefield based on the input-files of the system
	 */
	Spinorfield(const hardware::System&);

	/**
	 * Release resources
	 */
	~Spinorfield();

	/*
	 * Spinorfields cannot be copied
	 */
	Spinorfield& operator=(const Spinorfield&) = delete;
	Spinorfield(const Spinorfield&) = delete;
	Spinorfield() = delete;

	/**
	 * Get the buffers containing the gaugefield state on the devices.
	 */
	const std::vector<const hardware::buffers::Plain<spinor> *> get_buffers() const noexcept;

	/**
	 * Apply Gamma5 to the Spinorfield
	 */
	void gamma5() const;

	/**
	 * Set Spinorfield to zero
	 */
	void zero() const;

	/**
	 * Set Spinorfield to zero
	 */
	void cold() const;

	/**
	 * Set Spinorfield to be gaussian.
	 */
	void gaussian(const physics::PRNG& prng) const;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Plain<spinor> *> buffers;

	friend hmc_complex scalar_product(const Spinorfield& left, const Spinorfield& right);
	friend hmc_float squarenorm(const Spinorfield& field);
	friend void saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y);
	friend void sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x);
	friend void saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z);
};

/**
 * Create n spinorfields.
 *
 * \param n The number of spinorfields to create
 */
const std::vector<const Spinorfield *> create_spinorfields(const hardware::System& system, const size_t n)
;

/**
 * Release the given spinorfields
 */
void release_spinorfields(const std::vector<const Spinorfield *> fields);

/**
 * Calculate the scalar product of two spinorfields.
 */
hmc_complex scalar_product(const Spinorfield& left, const Spinorfield& right);
void scalar_product(const Scalar<hmc_complex>* res, const Spinorfield& left, const Spinorfield& right);

/**
 * Calculate the squarenorm of the spinorfield
 */
hmc_float squarenorm(const Spinorfield& field);
void squarenorm(const Scalar<hmc_float>* res, const Spinorfield& field);

/**
 * Perform the BLAS operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const Spinorfield& y);
void saxpy(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Spinorfield& y);

/**
 * Perform the BLAS operation sax.
 *
 * out = alpha * x
 */
void sax(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x);
void sax(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x);

/**
 * Perform the BLAS operation saxsbypz.
 *
 * out = alpha * x + beta * y + z
 */
void saxsbypz(const Spinorfield* out, const hmc_complex alpha, const Spinorfield& x, const hmc_complex beta, const Spinorfield& y, const Spinorfield& z);
void saxsbypz(const Spinorfield* out, const Scalar<hmc_complex>& alpha, const Spinorfield& x, const Scalar<hmc_complex>& beta, const Spinorfield& y, const Spinorfield& z);

}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_ */
