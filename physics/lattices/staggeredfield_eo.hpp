/** @file
 * Declaration of the physics::lattices::Staggeredfield_eo class
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

#ifndef _PHYSICS_LATTICES_STAGGEREDFIELD_EO_
#define _PHYSICS_LATTICES_STAGGEREDFIELD_EO_

#include "../../hardware/system.hpp"
#include "../../hardware/lattices/staggeredfield_eo.hpp"
#include "../../hardware/buffers/su3vec.hpp"
#include "../prng.hpp"
#include "scalar.hpp"
#include "vector.hpp"
#include "../../common_header_files/types_fermions.h"
#include "latticesInterfaces.hpp"
#include "../interfacesHandler.hpp"
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
	Staggeredfield_eo(const hardware::System&, const StaggeredfieldEoParametersInterface&);

	/**
	 * Release resources
	 */
	virtual ~Staggeredfield_eo(){}

	/**
	 * Staggeredfield_eo cannot be copied
	 */
	Staggeredfield_eo& operator=(const Staggeredfield_eo&) = delete;
	Staggeredfield_eo(const Staggeredfield_eo&) = delete;
	Staggeredfield_eo() = delete;
	Staggeredfield_eo(Staggeredfield_eo&&) = default;

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
	const StaggeredfieldEoParametersInterface& staggaredfieldEoParametersInterface;
	hardware::lattices::Staggeredfield_eo staggeredFieldEo;
	void import(const su3vec * const host) const;

	friend hmc_complex scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
	friend hmc_float scalar_product_real_part(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
	friend hmc_float squarenorm(const Staggeredfield_eo& field);
	friend void pseudo_randomize<Staggeredfield_eo, su3vec>(const Staggeredfield_eo* to, int seed);
};

/**
 * Create n Staggeredfield_eo
 */
std::vector<Staggeredfield_eo *> create_staggeredfields_eo(const hardware::System& system, const size_t n, physics::InterfacesHandler& interfacesHandler);

/**
 * Release the given Staggeredfield_eo
 */
void release_staggeredfields_eo(const std::vector<Staggeredfield_eo *> fields);

/**
 * Calculate the scalar product of two staggeredfields.
 */
hmc_complex scalar_product(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
void scalar_product(const Scalar<hmc_complex>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right);
hmc_float scalar_product_real_part(const Staggeredfield_eo& left, const Staggeredfield_eo& right);
void scalar_product_real_part(const Scalar<hmc_float>* res, const Staggeredfield_eo& left, const Staggeredfield_eo& right);

template<typename S, typename U, U (*T)(const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::scalar_product>(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::scalar_product_real_part>(const hardware::System&);

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
void sax(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x);
void sax(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x);
void sax(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x);

template<typename S, typename U, void (*T)(const S*, const U, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::sax>(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::sax>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
void saxpy(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
void saxpy(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
void saxpy(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);
void saxpy(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x, const Staggeredfield_eo& y);

template<typename S, typename U, void (*T)(const S*, const U, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpy>(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpy>(const hardware::System&);

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpby.
 *
 * out = alpha * x + beta * y
 */
void saxpby(const Staggeredfield_eo* out, const hmc_complex alpha, const Staggeredfield_eo& x, const hmc_complex beta, const Staggeredfield_eo& y);
void saxpby(const Staggeredfield_eo* out, const Scalar<hmc_complex>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_complex>& beta, const Staggeredfield_eo& y);
void saxpby(const Staggeredfield_eo* out, const hmc_float alpha, const Staggeredfield_eo& x, const hmc_float beta, const Staggeredfield_eo& y);
void saxpby(const Staggeredfield_eo* out, const Scalar<hmc_float>& alpha, const Staggeredfield_eo& x, const Scalar<hmc_float>& beta, const Staggeredfield_eo& y);
void saxpby(const Staggeredfield_eo* out, const Vector<hmc_float>& alpha, const int index_alpha, const Staggeredfield_eo& x, const Vector<hmc_float>& beta, const int index_beta, const Staggeredfield_eo& y);

template<typename S, typename U, void (*T)(const S*, const hmc_complex, const S&, const hmc_complex, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_complex, physics::lattices::saxpby>(const hardware::System&);
template<> size_t get_flops<physics::lattices::Staggeredfield_eo, hmc_float, physics::lattices::saxpby>(const hardware::System&);

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
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation sax + the squarenorm
 * with a set of real values of alpha.
 *
 * ||out||^2 = ||alpha[i] * x||^2
 */
void sax_vec_and_squarenorm(const Vector<hmc_float>* res, const Vector<hmc_float>& alpha, const Staggeredfield_eo& x);



/**
 * A utility function to log the tracenorm.
 *
 * It only evaluates in case the squarenorm will actually be printed.
 */
void log_squarenorm(const std::string& msg, const physics::lattices::Staggeredfield_eo& x);


}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_EO_ */

