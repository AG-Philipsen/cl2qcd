/** @file
 * physics::gaugeObservables class
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

#ifndef GAUGEOBSERVABLES_H_
#define GAUGEOBSERVABLES_H_

#include "../lattices/gaugefield.hpp"
#include "observablesInterfaces.hpp"

namespace physics{

  namespace observables{

    class Plaquettes
    {
    public:
      hmc_float plaquette;
      hmc_float temporalPlaquette;
      hmc_float spatialPlaquette;
      
    Plaquettes(hmc_float plaquetteIn, hmc_float temporalPlaquetteIn, hmc_float spatialPlaquetteIn):
      plaquette(plaquetteIn), temporalPlaquette(temporalPlaquetteIn), spatialPlaquette(spatialPlaquetteIn) {}
    };
    
    void measureGaugeObservablesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iteration, const physics::observables::GaugeObservablesParametersInterface&);
    hmc_float measurePlaquette(const physics::lattices::Gaugefield * gf, const physics::observables::GaugeObservablesParametersInterface&);
    hmc_float measurePlaquetteWithoutNormalization(const physics::lattices::Gaugefield * gf, const physics::observables::GaugeObservablesParametersInterface&);
    hmc_float measureRectangles(const physics::lattices::Gaugefield * gf, const physics::observables::GaugeObservablesParametersInterface&);
    Plaquettes measureAllPlaquettes(const physics::lattices::Gaugefield * gf, const physics::observables::GaugeObservablesParametersInterface&);
    hmc_complex measurePolyakovloop(const physics::lattices::Gaugefield * gf, const physics::observables::GaugeObservablesParametersInterface&);
  }
}

#endif

