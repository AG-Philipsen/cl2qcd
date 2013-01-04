/** @file
 * Declaration of the physics::lattices::Spinorfield class
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_
#define _PHYSICS_LATTICES_SPINORFIELD_

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

		private:
			hardware::System const& system;
			const std::vector<const hardware::buffers::Plain<spinor> *> buffers;
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

	}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_ */
