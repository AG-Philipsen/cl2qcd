/** @file
 * Declaration of the physics::lattices::Gaugemomenta class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_LATTICES_GAUGEMOMENTA_
#define _PHYSICS_LATTICES_GAUGEMOMENTA_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/gaugemomentum.hpp"
#include "../prng.hpp"
#include "scalar.hpp"

namespace physics {
namespace lattices {

template <class Lattice, typename Basetype> void pseudo_randomize(const Lattice* to, int seed);

/**
 * Representation of a gaugefield.
 */
class Gaugemomenta {

public:
	/**
	 * Construct a gaugemomenta field based on the input parameters of the system
	 */
	Gaugemomenta(const hardware::System&);

	/**
	 * Release resources
	 */
	virtual ~Gaugemomenta();

	/*
	 * Gaugemomenta cannot be copied
	 */
	Gaugemomenta& operator=(const Gaugemomenta&) = delete;
	Gaugemomenta(const Gaugemomenta&) = delete;
	Gaugemomenta() = delete;

	/**
	 * Get the buffers containing the gaugefield state on the devices.
	 */
	const std::vector<const hardware::buffers::Gaugemomentum *> get_buffers() const noexcept;

	/**
	 * Set the Gaugemomenta to zero
	 */
	void zero() const;

	/**
	 * Fill the Gaugemomenta according to a gaussian distribution.
	 *
	 * \exception std::invalid_argument if the prng is using different devices than this object
	 */
	void gaussian(const physics::PRNG& prng) const;

	/**
	 * Update the halo cells of each buffer from its neighbours.
	 *
	 * On a single device this will be a no-op.
	 */
	void update_halo() const;

	/**
	 * Get the number of elements.
	 */
	unsigned get_elements() const noexcept;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Gaugemomentum *> buffers;
	void import(const ae * const host) const;

	friend void saxpy(const Gaugemomenta* out, const hmc_float alpha, const Gaugemomenta& x, const Gaugemomenta& y);
	friend hmc_float squarenorm(const Gaugemomenta&);
	friend void pseudo_randomize<Gaugemomenta, ae>(const Gaugemomenta* to, int seed);
};

/**
 * Perform the BLAS (Basic Linear Algebra Subroutine) operation saxpy.
 *
 * out = alpha * x + y
 */
void saxpy(const Gaugemomenta* out, const hmc_float alpha, const Gaugemomenta& x, const Gaugemomenta& y);
void saxpy(const Gaugemomenta* out, const Scalar<hmc_float>& alpha, const Gaugemomenta& x, const Gaugemomenta& y);
//In the following two functions out = alpha * x + out
void saxpy(const Gaugemomenta* out, const hmc_float alpha, const Gaugemomenta& x);
void saxpy(const Gaugemomenta* out, const Scalar<hmc_float>& alpha, const Gaugemomenta& x);

template<typename S, void (*T)(const S*, const hmc_float, const S&, const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Gaugemomenta, physics::lattices::saxpy>(const hardware::System&);

/**
 * Calculate the squarenorm of the gaugemomenta
 */
hmc_float squarenorm(const Gaugemomenta& field);
void squarenorm(const Scalar<hmc_float>* res, const Gaugemomenta& field);

template<typename S, hmc_float (*T)(const S&)> size_t get_flops(const hardware::System&);
template<> size_t get_flops<physics::lattices::Gaugemomenta, physics::lattices::squarenorm>(const hardware::System&);

/**
 * A utility function to log the tracenorm.
 *
 * It only evaluates in case the squarenorm will actually be printed.
 */
void log_squarenorm(const std::string& msg, const physics::lattices::Gaugemomenta& x);
}
}

#endif /*_PHYSICS_LATTICES_GAUGEMOMENTA_ */
