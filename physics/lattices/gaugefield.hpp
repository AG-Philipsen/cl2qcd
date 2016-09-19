/** @file
 * Declaration of the physics::lattices::Gaugefield class
 *
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _PHYSICS_LATTICES_GAUGEFIELD_
#define _PHYSICS_LATTICES_GAUGEFIELD_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/su3.hpp"
#include "../prng.hpp"
#include "latticesInterfaces.hpp"
#include "../../hardware/lattices/gaugefield.hpp"

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
			Gaugefield(const hardware::System&, const GaugefieldParametersInterface * parameters, const physics::PRNG&);

			/**
			 * Construct a gaugefield based on the given ILDG file.
			 */
			Gaugefield(const hardware::System&, const GaugefieldParametersInterface * parameters, const physics::PRNG&, std::string);

			/**
			 * Construct a gaugefield that has been initialized hot or cold
			 */
			Gaugefield(const hardware::System&, const GaugefieldParametersInterface * parameters, const physics::PRNG&, bool hot);

			/**
			 * Release resources
			 */
			~Gaugefield();

			/*
			 * Gaugefields cannot be copied
			 */
			Gaugefield& operator=(const Gaugefield&) = delete;
			Gaugefield(const Gaugefield&) = delete;
			Gaugefield() = delete;

			/**
			 * Save gaugefield to a file with name conf_prefix + save + conf_postfix
			 * @param[in] number The trajectory number to be stored in the file
			 */
			void save(int number);
			/**
			 * Save gaugefield to a file with name conf_prefix + number + conf_postfix
			 * @param[in] number The trajectory number to be stored in the file
			 */
			void saveToSpecificFile(int number);
			/**
			 * Save gaugefield to a file with given name
			 * @param[in] outputfile name of file
			 * @param[in] number The trajectory number to be stored in the file
			 */
			void save(std::string outputfile, int number);

			/**
			 * Get the buffers containing the gaugefield state on the devices.
			 */
			const std::vector<const hardware::buffers::SU3 *> get_buffers() const noexcept;

			/**
			 * Smear the gaugefield.
			 * Creates a backup which can be restored via the unsmear command.
			 */
			void smear();
			void smear() const ;

			/**
			 * Restore the unsmeared gaugefield from the backup created before smearing.
			 */
			void unsmear();
			void unsmear() const;

			int get_trajectoryNumberAtInit() const;

			/**
			 * Update the halo cells of each buffer from its neighbours.
			 * On a single device this will be a no-op.
			 */
			void update_halo() const;

			std::string getName(int = -1) const noexcept;
			const physics::PRNG * getPrng() const;
			const hardware::System * getSystem() const;
			const GaugefieldParametersInterface * getParameters() const;
			
		private:
			hardware::System const& system;
			physics::PRNG const& prng;
			const GaugefieldParametersInterface * latticeObjectParameters;

			hardware::lattices::Gaugefield gaugefield;

			/**
			 * Utility functions for construction.
			 */
			void initializeBasedOnParameters();
			void initializeHotOrCold(bool hot);
			void initializeFromILDGSourcefile(std::string);

			int trajectoryNumberAtInit;
		};

	}
}



#endif /*_PHYSICS_LATTICES_GAUGEFIELD_ */
