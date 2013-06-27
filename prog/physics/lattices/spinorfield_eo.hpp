/** @file
 * Declaration of the physics::lattices::Spinorfield_eo class
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_EO_
#define _PHYSICS_LATTICES_SPINORFIELD_EO_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/spinor.hpp"
#include "../prng.hpp"
#include "spinorfield.hpp"
#include "scalar.hpp"
#include "../../types_fermions.h"

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

	/**
	 * Mark the halo as requiring an update.
	 *
	 * @warning: The transfer may be defferred. Make sure to call require_halo if you require any halo values.
	 */
	void mark_halo_dirty() const;

	/**
	 * Ensure that the halo has been update.
	 *
	 * \param width Only require the the given width of the halo to be up to date, use 0 to indicate the full halo is required.
	 */
	void require_halo(unsigned width = 0) const;

	/**
	 * Mark the halo as up to date.
	 */
	void mark_halo_clean() const;


private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Spinor *> buffers;

	void update_halo() const;
#ifdef LAZY_HALO_UPDATES
	mutable bool halo_dirty;
#endif

	friend hmc_complex scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right);
	friend hmc_float squarenorm(const Spinorfield_eo& field);
	friend void saxpy(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const Spinorfield_eo& y);
	friend void sax(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x);
	friend void saxsbypz(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const hmc_complex beta, const Spinorfield_eo& y, const Spinorfield_eo& z);
};

/**
 * Copy the contents of one lattice to another
 *
 * \param[out] dest The lattice to copy to
 * \param[in]  from The lattice to copy from
 */
void copyData(const Spinorfield_eo* to, const Spinorfield_eo* from);

/**
 * Copy the contents of one lattice to another
 *
 * \param[out] dest Spinorfield_eohe lattice to copy to
 * \param[in]  from Spinorfield_eohe lattice to copy from
 */
void copyData(const Spinorfield_eo* to, const Spinorfield_eo& from);

/**
 * Calculate the scalar product of two spinorfields.
 */
hmc_complex scalar_product(const Spinorfield_eo& left, const Spinorfield_eo& right);
void scalar_product(const Scalar<hmc_complex>* res, const Spinorfield_eo& left, const Spinorfield_eo& right);

template<typename S, hmc_complex (*T)(const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Spinorfield_eo, physics::lattices::scalar_product>(const hardware::System&);

/**
 * Calculate the squarenorm of the spinorfield
 */
hmc_float squarenorm(const Spinorfield_eo& field);
void squarenorm(const Scalar<hmc_float>* res, const Spinorfield_eo& field);

template<typename S, hmc_float (*T)(const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Spinorfield_eo, physics::lattices::squarenorm>(const hardware::System&);

/**
 * Perform the BLAS operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const Spinorfield_eo& y);
void saxpy(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x, const Spinorfield_eo& y);

template<typename S, void (*T)(const S*, const hmc_complex, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxpy>(const hardware::System&);

/**
 * Perform the BLAS operation sax.
 *
 * out = alpha * x
 */
void sax(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x);
void sax(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x);

template<typename S, void (*T)(const S*, const hmc_complex, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Spinorfield_eo, physics::lattices::sax>(const hardware::System&);

/**
 * Perform the BLAS operation saxsbypz.
 *
 * out = alpha * x + beta * y + z
 */
void saxsbypz(const Spinorfield_eo* out, const hmc_complex alpha, const Spinorfield_eo& x, const hmc_complex beta, const Spinorfield_eo& y, const Spinorfield_eo& z);
void saxsbypz(const Spinorfield_eo* out, const Scalar<hmc_complex>& alpha, const Spinorfield_eo& x, const Scalar<hmc_complex>& beta, const Spinorfield_eo& y, const Spinorfield_eo& z);

template<typename S, void (*T)(const S*, const hmc_complex, const S&, const hmc_complex, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Spinorfield_eo, physics::lattices::saxsbypz>(const hardware::System&);

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

/**
 * A utility function to log the tracenorm.
 *
 * It only evaluates in case the squarenorm will actually be printed.
 */
void log_squarenorm(const std::string& msg, const physics::lattices::Spinorfield_eo& x);
}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_EO_ */
