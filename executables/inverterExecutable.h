/** @file
 *
 * Everything required by inverter's main()
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
#ifndef _INVERTERH_
#define _INVERTERH_

#include "measurementExecutable.h"
#include "../physics/lattices/spinorfield.hpp"
#include "../physics/sources.hpp"
#include "../physics/algorithms/flavour_doublet.hpp"
#include "../physics/algorithms/inversion.hpp"

/**
 * Inverter executable which measures fermionic observables on gaugefield configurations.
 */
class inverterExecutable : public measurementExecutable
{
public:
	/**
	 * Constructor.
	 * Writes Logfile.
	 */
	inverterExecutable(int argc, const char* argv[]);

protected:
	const std::string 	filenameForInverterLogfile = "inverter.log";
	std::string 		filenameForTwoFlavourDoubletChiralCondensateData;
	std::string 		filenameForTwoFlavourDoubletCorrelatorData;

	void writeInverterLogfile();

	void printParametersToScreenAndFile();

	void measureTwoFlavourDoubletCorrelatorsOnGaugefield();

	void measureTwoFlavourDoubletChiralCondensateOnGaugefield();

	/**
	 * Performs measurements of fermionic observables on possibly multiple gaugefield configurations.
	 */
	void performApplicationSpecificMeasurements();
};

#endif /* _INVERTERH_ */

