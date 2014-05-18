/** @file
 * physics::gaugeObservables class
 *
 * Copyright 2014,Christopher Pinke
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

  class gaugeObservables{
  public:
    gaugeObservables(meta::Inputparameters * parametersIn)
      {
	parameters = parametersIn;
      }
    gaugeObservables() = delete;

    /**
     * Measures all gauge observables according to parameter settings
     */
    void measureGaugeObservables(physics::lattices::Gaugefield& gf, int iteration, meta::Inputparameters params);

    /**
     * Measures plaquette and polyakov loop and writes them to file
     */
    void measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params);
    void measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params);

    /**
     * Measures rectangle products of link variables and writes it to file
     */
    void measureRectangles(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params);

    /**
     * Measures transportcoefficient kappa and writes it to file
     */
    void measureTransportcoefficientKappa(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters);

  private:
    meta::Inputparameters * parameters;
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
    double getPlaquette()
    {
      return plaquette;
    }

    void measurePlaquette(physics::lattices::Gaugefield * gaugefield)
    {
      //      void physics::lattices::Gaugefield::plaquette(hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq) const
	// the plaquette is local to each side and then summed up
	// for multi-device simply calculate the plaquette for each device and then sum up the devices

	using hardware::buffers::Plain;

	auto gaugefieldBuffers = gaugefield->get_buffers();
	size_t num_devs = gaugefieldBuffers.size();

	if(num_devs == 1) {
		auto device = gaugefieldBuffers[0]->get_device();

		const Plain<hmc_float> plaq_dev(1, device);
		const Plain<hmc_float> tplaq_dev(1, device);
		const Plain<hmc_float> splaq_dev(1, device);
		device->get_gaugefield_code()->plaquette_device(gaugefieldBuffers[0], &plaq_dev, &tplaq_dev, &splaq_dev);

		plaq_dev.dump(&plaquette);
		tplaq_dev.dump(&plaquette_temporal);
		splaq_dev.dump(&plaquette_spatial);
		plaquette_temporal /= static_cast<hmc_float>(meta::get_tplaq_norm(*parameters));
		plaquette_spatial /= static_cast<hmc_float>(meta::get_splaq_norm(*parameters));
		plaquette  /= static_cast<hmc_float>(meta::get_plaq_norm(*parameters));
	} else {
		// trigger calculation
		std::vector<const Plain<hmc_float>*> plaqs; plaqs.reserve(num_devs);
		std::vector<const Plain<hmc_float>*> tplaqs; tplaqs.reserve(num_devs);
		std::vector<const Plain<hmc_float>*> splaqs; splaqs.reserve(num_devs);
		for(size_t i = 0; i < num_devs; ++i) {
			auto device = gaugefieldBuffers[i]->get_device();
			const Plain<hmc_float>* plaq_dev = new Plain<hmc_float>(1, device);
			const Plain<hmc_float>* tplaq_dev = new Plain<hmc_float>(1, device);
			const Plain<hmc_float>* splaq_dev = new Plain<hmc_float>(1, device);
			device->get_gaugefield_code()->plaquette_device(gaugefieldBuffers[i], plaq_dev, tplaq_dev, splaq_dev);
			plaqs.push_back(plaq_dev);
			tplaqs.push_back(tplaq_dev);
			splaqs.push_back(splaq_dev);
		}
		// collect results
		plaquette = 0.0;
		plaquette_spatial = 0.0;
		plaquette_temporal = 0.0;
		for(size_t i = 0; i < num_devs; ++i) {
			double tmp;

			plaqs[i]->dump(&tmp);
			plaquette += tmp;

			tplaqs[i]->dump(&tmp);
			plaquette_temporal += tmp;

			splaqs[i]->dump(&tmp);
			plaquette_spatial += tmp;

			delete plaqs[i];
			delete tplaqs[i];
			delete splaqs[i];
		}
		plaquette_temporal /= static_cast<hmc_float>(meta::get_tplaq_norm(*parameters));
		plaquette_spatial /= static_cast<hmc_float>(meta::get_splaq_norm(*parameters));
		plaquette  /= static_cast<hmc_float>(meta::get_plaq_norm(*parameters));
	}
    }
  };
}

#endif
