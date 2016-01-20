/** @file
 * Implemenation of physics::gaugeObservables class
 *
 * Copyright 2013 Matthias Bach
 *           2014 Christopher Pinke
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

#include "gaugeObservables.hpp"
#include "../../hardware/code/kappa.hpp"
#include <cassert>
#include <fstream>
#include <cmath>
#include "../../hardware/code/gaugefield.hpp"

class gaugeObservables{
    public:
        gaugeObservables(const physics::observables::GaugeObservablesParametersInterface& interface)
            : gaugeObservablesParametersInterface(interface), outputToFile("") {};
        gaugeObservables() = delete;
        void measureGaugeObservablesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iteration);
        void measurePlaqAndPolyAndWriteToFile(const physics::lattices::Gaugefield * gf, int iter, const std::string& filename);
        void measurePlaqAndPolyAndWriteToFile(const physics::lattices::Gaugefield * gf, int iter);
        void measureRectanglesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iter);
        void measureTransportcoefficientKappaAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, int iteration);
        physics::observables::Plaquettes measurePlaquettes(const physics::lattices::Gaugefield * gaugefield, bool normalize = true);
        hmc_float measureRectangles(const physics::lattices::Gaugefield * gaugefield);
        hmc_complex measurePolyakovloop(const physics::lattices::Gaugefield * gaugefield);

    private:
        const physics::observables::GaugeObservablesParametersInterface& gaugeObservablesParametersInterface;
        std::ofstream outputToFile;
        void writePlaqAndPolyToFile(int iter,  const std::string& filename, const physics::observables::Plaquettes plaquettes, const hmc_complex polyakov);
        void writeTransportcoefficientKappaToFile(std::string filename, int iteration, const hmc_float kappa);
        void writeTransportcoefficientKappaToFileUsingOpenOutputStream(int iteration, const hmc_float kappa);
        void writeRectanglesToFile(int iter, const std::string& filename, hmc_float rectangles);
};

void gaugeObservables::measureGaugeObservablesAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, int iteration)
{
  measurePlaqAndPolyAndWriteToFile(gaugefield, iteration);
  if ( gaugeObservablesParametersInterface.measureRectangles() ){
    measureRectanglesAndWriteToFile(gaugefield, iteration);
  }
  if ( gaugeObservablesParametersInterface.measureTransportCoefficientKappa() ) {
    measureTransportcoefficientKappaAndWriteToFile(gaugefield, iteration);
  }
}

void gaugeObservables::measurePlaqAndPolyAndWriteToFile(const physics::lattices::Gaugefield * gf, int iter)
{
	const std::string filename = gaugeObservablesParametersInterface.getGaugeObservablesFilename("");
	measurePlaqAndPolyAndWriteToFile(gf, iter, filename);
}

void gaugeObservables::measurePlaqAndPolyAndWriteToFile(const physics::lattices::Gaugefield * gf, int iter, const std::string& filename)
{
  physics::observables::Plaquettes plaquettes = measurePlaquettes(gf);
  hmc_complex polyakov = measurePolyakovloop(gf);
	if ( gaugeObservablesParametersInterface.printToScreen() ){
	  logger.info() << iter << '\t' << plaquettes.plaquette << '\t' << plaquettes.temporalPlaquette << '\t' << plaquettes.spatialPlaquette << '\t' << polyakov.re << '\t' << polyakov.im << '\t' << sqrt(polyakov.re * polyakov.re + polyakov.im * polyakov.im);
	}
	writePlaqAndPolyToFile(iter, filename, plaquettes, polyakov);
}

void gaugeObservables::writePlaqAndPolyToFile(int iter,  const std::string& filename, const physics::observables::Plaquettes plaquettes, const hmc_complex polyakov)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile.width(8);
	outputToFile << iter << "\t";
	outputToFile.precision(15);
	outputToFile << plaquettes.plaquette << "\t" << plaquettes.temporalPlaquette<< "\t" << plaquettes.spatialPlaquette<< "\t" << polyakov.re << "\t" << polyakov.im << "\t" << sqrt(polyakov.re * polyakov.re + polyakov.im * polyakov.im) << std::endl;
	outputToFile.close();
}

static hmc_float kappa_clover(const physics::lattices::Gaugefield& gf, hmc_float beta)
{
	assert(gf.get_buffers().size() == 1);

	auto gf_dev = gf.get_buffers()[0];

	auto device = gf_dev->get_device();
	hardware::buffers::Plain<hmc_float> kappa_clover_dev(1, device);
	gf_dev->get_device()->getKappaCode()->run_kappa_clover(&kappa_clover_dev, gf_dev, beta);

	hmc_float kappa_clover_host;
	kappa_clover_dev.dump(&kappa_clover_host);
	return kappa_clover_host;
}

void gaugeObservables::measureTransportcoefficientKappaAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, int iteration)
{
	hmc_float kappa = kappa_clover(*gaugefield, gaugeObservablesParametersInterface.getBeta());
	writeTransportcoefficientKappaToFile(gaugeObservablesParametersInterface.getTransportCoefficientKappaFilename(), iteration, kappa);
}

void gaugeObservables::writeTransportcoefficientKappaToFileUsingOpenOutputStream(int iteration, const hmc_float kappa)
{
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << iteration << "\t" << kappa << std::endl;
}

void gaugeObservables::writeTransportcoefficientKappaToFile(std::string filename, int iteration, const hmc_float kappa)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if ( outputToFile.is_open() ) {
	  writeTransportcoefficientKappaToFileUsingOpenOutputStream(iteration, kappa);
	  outputToFile.close();
	} else {
		logger.warn() << "Could not open " << filename;
		File_Exception(filename.c_str());
	}
}

void gaugeObservables::measureRectanglesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iteration)
{
	hmc_float rectangles = measureRectangles(gf);
	writeRectanglesToFile(iteration, gaugeObservablesParametersInterface.getRectanglesFilename(), rectangles);
}

void gaugeObservables::writeRectanglesToFile(int iter, const std::string& filename, const hmc_float rectangles)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile.width(8);
	outputToFile << iter << "\t";
	outputToFile.precision(15);
	outputToFile << rectangles << std::endl;
	outputToFile.close();
}

physics::observables::Plaquettes gaugeObservables::measurePlaquettes(const physics::lattices::Gaugefield * gaugefield, bool normalize)
{
    // the plaquette is local to each side and then summed up
    // for multi-device simply calculate the plaquette for each device and then sum up the devices

    using hardware::buffers::Plain;

    physics::observables::Plaquettes plaquettes{0.0, 0.0, 0.0};
    auto gaugefieldBuffers = gaugefield->get_buffers();
    size_t num_devs = gaugefieldBuffers.size();

    if(num_devs == 1) {
        auto device = gaugefieldBuffers[0]->get_device();

        const Plain<hmc_float> plaq_dev(1, device);
        const Plain<hmc_float> tplaq_dev(1, device);
        const Plain<hmc_float> splaq_dev(1, device);
		device->getGaugefieldCode()->plaquette_device(gaugefieldBuffers[0], &plaq_dev, &tplaq_dev, &splaq_dev);

        plaq_dev.dump(&plaquettes.plaquette);
        tplaq_dev.dump(&plaquettes.temporalPlaquette);
        splaq_dev.dump(&plaquettes.spatialPlaquette);
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
			device->getGaugefieldCode()->plaquette_device(gaugefieldBuffers[i], plaq_dev, tplaq_dev, splaq_dev);
            plaqs.push_back(plaq_dev);
            tplaqs.push_back(tplaq_dev);
            splaqs.push_back(splaq_dev);
        }
        // collect results
        plaquettes.plaquette = 0.0;
        plaquettes.spatialPlaquette = 0.0;
        plaquettes.temporalPlaquette = 0.0;
        for(size_t i = 0; i < num_devs; ++i) {
            hmc_float tmp;

            plaqs[i]->dump(&tmp);
            plaquettes.plaquette += tmp;

            tplaqs[i]->dump(&tmp);
            plaquettes.temporalPlaquette += tmp;

            splaqs[i]->dump(&tmp);
            plaquettes.spatialPlaquette += tmp;

            delete plaqs[i];
            delete tplaqs[i];
            delete splaqs[i];
        }
    }

    if (normalize)
    {
        plaquettes.temporalPlaquette/= static_cast<hmc_float>(gaugeObservablesParametersInterface.getTemporalPlaquetteNormalization());
        plaquettes.spatialPlaquette /= static_cast<hmc_float>(gaugeObservablesParametersInterface.getSpatialPlaquetteNormalization());
        plaquettes.plaquette  /= static_cast<hmc_float>(gaugeObservablesParametersInterface.getPlaquetteNormalization());
    }

    return plaquettes;
}

hmc_float gaugeObservables::measureRectangles(const physics::lattices::Gaugefield * gaugefield)
{
    // the rectangles are local to each site and then summed up
    // for multi-device simply calculate the plaquette for each device and then sum up the devices

    using hardware::buffers::Plain;

    hmc_float rectangles;
    auto gaugefieldBuffers = gaugefield->get_buffers();
    size_t num_devs = gaugefieldBuffers.size();

    if(num_devs == 1) {
        auto device = gaugefieldBuffers[0]->get_device();
        const Plain<hmc_float> rect_dev(1, device);
        device->getGaugefieldCode()->rectangles_device(gaugefieldBuffers[0], &rect_dev);
        rect_dev.dump(&rectangles);
    } else {
        // trigger calculation
        std::vector<const Plain<hmc_float>*> rects;
        rects.reserve(num_devs);
        for(size_t i = 0; i < num_devs; ++i) {
            auto device = gaugefieldBuffers[i]->get_device();
            const Plain<hmc_float>* rect_dev = new Plain<hmc_float>(1, device);
            device->getGaugefieldCode()->rectangles_device(gaugefieldBuffers[i], rect_dev);
            rects.push_back(rect_dev);
        }
        // collect results
        rectangles = 0.0;
        for(size_t i = 0; i < num_devs; ++i) {
            hmc_float tmp;
            rects[i]->dump(&tmp);
            rectangles += tmp;

            delete rects[i];
        }
    }
    return rectangles;
}

hmc_complex gaugeObservables::measurePolyakovloop(const physics::lattices::Gaugefield * gaugefield)
{
    using hardware::buffers::Plain;

    hmc_complex polyakov;
    auto gaugefieldBuffers = gaugefield->get_buffers();
    size_t num_devs = gaugefieldBuffers.size();

    if(num_devs == 1) {
        auto gf_buf = gaugefieldBuffers[0];
        auto device = gf_buf->get_device();
        auto gf_code = device->getGaugefieldCode();
        const Plain<hmc_complex> pol_buf(1, device);

        gf_code->polyakov_device(gf_buf, &pol_buf);
        pol_buf.dump(&polyakov);
    } else {
        const size_t volspace = gaugeObservablesParametersInterface.getSpatialVolume();
        // calculate local part per device
        std::vector<Plain<Matrixsu3>*> local_results;
        local_results.reserve(num_devs);
        for(auto buffer: gaugefieldBuffers) {
            auto dev = buffer->get_device();
            auto res_buf = new Plain<Matrixsu3>(volspace, dev);
            dev->getGaugefieldCode()->polyakov_md_local_device(res_buf, buffer);
            local_results.push_back(res_buf);
        }

        // merge results
        auto main_dev = gaugefieldBuffers.at(0)->get_device();
        Plain<Matrixsu3> merged_buf(volspace * local_results.size(), main_dev);
        {
            Matrixsu3* merged_host = new Matrixsu3[merged_buf.get_elements()];
            size_t offset = 0;
            for(auto local_res: local_results) {
                local_res->dump(&merged_host[offset]);
                offset += volspace;
            }
            merged_buf.load(merged_host);
            delete[] merged_host;
        }

        const Plain<hmc_complex> pol_buf(1, main_dev);
        main_dev->getGaugefieldCode()->polyakov_md_merge_device(&merged_buf, num_devs, &pol_buf);
        pol_buf.dump(&polyakov);

        for(auto buffer: local_results) {
            delete buffer;
        }
    }

    polyakov.re /= static_cast<hmc_float>(gaugeObservablesParametersInterface.getPolyakovLoopNormalization());
    polyakov.im /= static_cast<hmc_float>(gaugeObservablesParametersInterface.getPolyakovLoopNormalization());

    return polyakov;
}

void physics::observables::measureGaugeObservablesAndWriteToFile(const physics::lattices::Gaugefield * gf, int iteration,
                                                                 const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    obs.measureGaugeObservablesAndWriteToFile(gf, iteration);
}

hmc_float physics::observables::measurePlaquette(const physics::lattices::Gaugefield * gf,
                                                 const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    return obs.measurePlaquettes(gf).plaquette;
}

hmc_float physics::observables::measurePlaquetteWithoutNormalization(const physics::lattices::Gaugefield * gf,
                                                                     const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    return obs.measurePlaquettes(gf, false).plaquette;
}

hmc_float physics::observables::measureRectangles(const physics::lattices::Gaugefield * gf,
                                                  const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    return obs.measureRectangles(gf);
}

hmc_complex  physics::observables::measurePolyakovloop(const physics::lattices::Gaugefield * gf,
                                                       const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    return obs.measurePolyakovloop(gf);
}

physics::observables::Plaquettes  physics::observables::measureAllPlaquettes(const physics::lattices::Gaugefield * gf,
                                                                             const physics::observables::GaugeObservablesParametersInterface& parameters)
{
    gaugeObservables obs(parameters);
    return obs.measurePlaquettes(gf);
}
