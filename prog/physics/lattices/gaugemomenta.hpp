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
	 * Gaugemomentas cannot be copied
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

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Gaugemomentum *> buffers;

	friend hmc_float squarenorm(const Gaugemomenta&);
};

/**
 * Calculate the squarenorm of the gaugemomenta
 */
hmc_float squarenorm(const Gaugemomenta& field);
void squarenorm(const Scalar<hmc_float>* res, const Gaugemomenta& field);

}
}

#endif /*_PHYSICS_LATTICES_GAUGEMOMENTA_ */
