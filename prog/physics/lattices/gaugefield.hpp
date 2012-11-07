/** @file
 * Declaration of the physics::lattices::Gaugefield class
 */

#ifndef _PHYSICS_LATTICES_GAUGEFIELD_
#define _PHYSICS_LATTICES_GAUGEFIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/su3.hpp"
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
		class Gaugefield {

		public:
			/**
			 * Construct a gaugefield based on the input-files of the system
			 */
			Gaugefield(hardware::System&, physics::PRNG&);

			/**
			 * Construct a gaugefield based on the given ILDG file.
			 */
			Gaugefield(hardware::System&, physics::PRNG&, std::string);

			/**
			 * Construct a gaugefield that has been initialized hot or cold
			 */
			Gaugefield(hardware::System&, physics::PRNG&, bool hot);

			/*
			 * Gaugefields cannot be copied
			 */
			Gaugefield& operator=(const Gaugefield&) = delete;
			Gaugefield(const Gaugefield&) = delete;
			Gaugefield() = delete;

			/**
			 * Save gaugefield to a file with name conf.number
			 * @param[in] number number to be added to file name
			 */
			void save(int number);
			/**
			 * Save gaugefield to a file with given name
			 * @param[in] outputfile name of file
			 */
			void save(std::string outputfile);

			/**
			 * Return plaquette value
			 */
			hmc_float plaquette();

		private:
			hardware::System const& system;
			physics::PRNG const& prng;
			const std::vector<const hardware::buffers::SU3 *> buffers;

			/**
			 * Utility function for construction.
			 */
			void fill_from_ildg(std::string);
		};
	}
}

#endif /*_PHYSICS_LATTICES_GAUGEFIELD_ */
