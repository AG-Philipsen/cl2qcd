/** @file
 * Declaration of the physics::lattices::Spinorfield class
 */

#ifndef _PHYSICS_LATTICES_GAUGEFIELD_
#define _PHYSICS_LATTICES_GAUGEFIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/plain.hpp"
#include "../prng.hpp"

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
			Spinorfield(hardware::System&);

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

		private:
			hardware::System const& system;
			const std::vector<const hardware::buffers::Plain<spinor> *> buffers;
		};
	}
}

#endif /*_PHYSICS_LATTICES_GAUGEFIELD_ */
