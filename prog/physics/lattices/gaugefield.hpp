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
			hmc_float plaquette() const;

			/**
			 * Calculate plaquette and polyakov of this gaugefield.
			 *
			 * @param[out] plaq Storage for result of plaquette calculation
			 * @param[out] tplaq Storage for result of plaquette calculation
			 * @param[out] splaq Storage for result of plaquette calculation
			 * @param[out] pol Storage for result of polyakov calculation
			 */
			void gaugeobservables(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol) const;

			/**
			 * Calculate rectangles of this gaugefield
			 *
			 * @param[in] gf gaugefield to measure on
			 * @param[out] plaq Storage for result of rectangles calculation
			 */
			hmc_float rectangles() const;

			/**
			 * Get the buffers containing the gaugefield state on the devices.
			 */
			const std::vector<const hardware::buffers::SU3 *> get_buffers() const noexcept;

		private:
			hardware::System const& system;
			physics::PRNG const& prng;
			const std::vector<const hardware::buffers::SU3 *> buffers;

			/**
			 * Utility function for construction.
			 */
			void fill_from_ildg(std::string);
		};

		/**
		 * Print the gaugeobservables of the given gaugefield as info output.
		 */
		void print_gaugeobservables(const physics::lattices::Gaugefield& gf, int iter);

		/**
		 * Print the gaugeobservables of the given gaugefield to the given file
		 */
		void print_gaugeobservables(const physics::lattices::Gaugefield& gf, int iter, const std::string& filename);
	}
}

#endif /*_PHYSICS_LATTICES_GAUGEFIELD_ */
