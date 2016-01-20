/** @file
 * physics::observables::wilson::TwoFlavourChiralCondensate class
 *
 * Copyright 2014 Christopher Pinke
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

#ifndef WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_
#define WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include <fstream>
#include "observablesInterfaces.hpp"

namespace physics {

    namespace observables {

        namespace wilson {

            std::vector<double> measureTwoFlavourChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, std::string configurationName,
                                                                                physics::InterfacesHandler & interfacesHandler);
            std::vector<double> measureTwoFlavourChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, int iteration,
                                                                                physics::InterfacesHandler & interfacesHandler);
        }
    }
}
#endif
