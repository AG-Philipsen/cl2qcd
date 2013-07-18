/** @file
 * Declaration of the physics::lattices::Staggeredfield_eo class
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#ifndef _PHYSICS_LATTICES_STAGGEREDFIELD_EO_
#define _PHYSICS_LATTICES_STAGGEREDFIELD_EO_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/su3vec.hpp"
#include "../prng.hpp"
//#include "spinorfield.hpp"
#include "scalar.hpp"
#include "../../types_fermions.h"
//This is to make the template pseudo_randomize friend of this class
#include "util.hpp"


/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a staggeredfield (with eo preconditioning).
 */
class Staggeredfield_eo {

public:
	/**
	 * Construct a staggeredfield based on the input-files of the system
	 */
	Staggeredfield_eo(const hardware::System&);

	/**
	 * Release resources
	 */
	~Staggeredfield_eo();

	/**
	 * Staggeredfield_eo cannot be copied
	 */
	Staggeredfield_eo& operator=(const Staggeredfield_eo&) = delete;
	Staggeredfield_eo(const Staggeredfield_eo&) = delete;
	Staggeredfield_eo() = delete;

	/**
	 * Get the buffers containing the staggeredfield state on the devices.
	 */
	const std::vector<const hardware::buffers::SU3vec *> get_buffers() const noexcept;

	/**
	 * Set Staggeredfield to zero
	 */
	void set_zero() const;

	/**
	 * Set Staggeredfield to cold
	 */
	void set_cold() const;

	/**
	 * Set Staggeredfield to be gaussian.
	 */
	void set_gaussian(const physics::PRNG& prng) const;

	/**
	 * Update the halos of the spinorfield buffers.
	 */
	void update_halo() const;
	
	/**
	 * Get the number of elements of the field, namely the lattice volume.
	 */
	unsigned get_elements() const noexcept;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::SU3vec *> buffers;
	void import(const su3vec * const host) const;

	friend hmc_complex scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
	friend hmc_float squarenorm(const Staggeredfield_eo& field);
	friend void sax(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x);
	friend void saxpy(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
	friend void saxpby(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y);
	friend void saxpbypz(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z);
	friend void pseudo_randomize<Staggeredfield_eo, su3vec>(const Staggeredfield_eo* to, int seed);
};


/**
 * Calculate the scalar product of two staggeredfields.
 */
hmc_complex scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
void scalar_product(const Scalar<hmc_complex>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right);

template<typename S, hmc_complex (*T)(const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::scalar_product>(const hardware::System&);

/**
 * Calculate the squarenorm of the spinorfield
 */
hmc_float squarenorm(const Staggeredfield_eo& field);
void squarenorm(const Scalar<hmc_float>* res, const Staggeredfield_eo& field);

template<typename S, hmc_float (*T)(const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::squarenorm>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation sax.
 *
 * out = alpha * x
 */
void sax(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x);
void sax(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x);

template<typename S, void (*T)(const S*, const hmc_complex, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::sax>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
void saxpy(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);

template<typename S, void (*T)(const S*, const hmc_complex, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::saxpy>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpby.
 *
 * out = alpha * x + beta * y
 */
void saxpby(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y);
void saxpby(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y);

template<typename S, void (*T)(const S*, const hmc_complex, const S&, const hmc_complex, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::saxpby>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpbypz.
 *
 * out = alpha * x + beta * y + z
 */
void saxpbypz(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z);
void saxpbypz(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y, const Staggeredfield_eo& z);

template<typename S, void (*T)(const S*, const hmc_complex, const S&, const hmc_complex, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, physics::lattices::saxpbypz>(const hardware::System&);

/**
 * A utility function to log the tracenorm.
 *
 * It only evaluates in case the squarenorm will actually be printed.
 */
void log_squarenorm(const std::string& msg, const physics::lattices::Staggeredfield_eo& x);



//So far not implemented, since only EO-preconditioning is used
#if 0
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
#endif

}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_EO_ */

