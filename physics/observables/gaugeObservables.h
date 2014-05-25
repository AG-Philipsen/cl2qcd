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
#include <fstream>
#include <cmath>
#include "../../meta/inputparameters.hpp"
#include "../algorithms/kappa_clover.hpp"

namespace physics{

  namespace observables{

    class Plaquettes
    {
    public:
      double plaquette;
      double temporalPlaquette;
      double spatialPlaquette;
      
    Plaquettes(double plaquetteIn, double temporalPlaquetteIn, double spatialPlaquetteIn):
      plaquette(plaquetteIn), temporalPlaquette(temporalPlaquetteIn), spatialPlaquette(spatialPlaquetteIn) {}
    };
    
    class gaugeObservables{
    public:
      gaugeObservables(const meta::Inputparameters * parametersIn)
	{
	  parameters = parametersIn;
	}
      gaugeObservables() = delete;
      
      /**
       * Measures all gauge observables according to parameter settings
       */
      void measureGaugeObservables(const physics::lattices::Gaugefield * gf, int iteration);
      
      /**
       * Measures plaquette and polyakov loop and writes them to file
       */
      void measurePlaqAndPoly(const physics::lattices::Gaugefield * gf, int iter, const std::string& filename);
      void measurePlaqAndPoly(const physics::lattices::Gaugefield * gf, int iter);
      
      /**
       * Measures rectangle products of link variables and writes it to file
       */
      void measureRectangles(const physics::lattices::Gaugefield * gf, int iter);
      
      /**
       * Measures transportcoefficient kappa and writes it to file
       */
      void measureTransportcoefficientKappa(const physics::lattices::Gaugefield * gaugefield, int iteration);
      
    private:
      const meta::Inputparameters * parameters;
      double plaquette;
      double plaquette_temporal;
      double plaquette_spatial;
      double rectangles;
      double kappa;
      hmc_complex polyakov;
      std::ofstream outputToFile;
      
      void writePlaqAndPolyToFile(int iter,  const std::string& filename);
      
      void writeTransportcoefficientKappaToFile(std::string filename, int iteration);
      
      void writeTransportcoefficientKappaToFileUsingOpenOutputStream(int iteration);
      
      void writeRectanglesToFile(int iter, const std::string& filename);
      
    public:
      double getPlaquette();
      double getSpatialPlaquette();
      double getTemporalPlaquette();
      double getRectangles();
      hmc_complex getPolyakovloop();
      
      double measurePlaquette(const physics::lattices::Gaugefield * gaugefield);
      double measureRectangles(const physics::lattices::Gaugefield * gaugefield);
      hmc_complex measurePolyakovloop(const physics::lattices::Gaugefield * gaugefield);
    };

    void measureGaugeObservablesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iteration);
    double measurePlaquette(const physics::lattices::Gaugefield * gf);
    double measureRectangles(const physics::lattices::Gaugefield * gf);
    Plaquettes measureAllPlaquettes(const physics::lattices::Gaugefield * gf);
    hmc_complex measurePolyakovloop(const physics::lattices::Gaugefield * gf);
  }
}

#endif
