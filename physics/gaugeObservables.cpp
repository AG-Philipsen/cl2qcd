#include "gaugeObservables.h"

void physics::gaugeObservables::measureGaugeObservables(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters)
{
  measurePlaqAndPoly(gaugefield, iteration, parameters);
  if ( parameters.get_measure_rectangles() ){
    measureRectangles(gaugefield, iteration, parameters);
  }
  if ( parameters.get_measure_transportcoefficient_kappa() ) {
    measureTransportcoefficientKappa(gaugefield, iteration, parameters);
  }
}

void physics::gaugeObservables::measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params)
{
	const std::string filename = meta::get_gauge_obs_file_name(params,  "");
	measurePlaqAndPoly(gf, iter, filename, params);
}

void physics::gaugeObservables::measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params)
{
	gf.gaugeobservables(&plaquette, &plaquette_temporal, &plaquette_spatial, &polyakov);

	if ( params.get_print_to_screen() ){
	  logger.info() << iter << '\t' << plaquette << '\t' << plaquette_temporal << '\t' << plaquette_spatial << '\t' << polyakov.re << '\t' << polyakov.im << '\t' << sqrt(polyakov.re * polyakov.re + polyakov.im * polyakov.im);
	}
	writePlaqAndPolyToFile(iter, filename);
}

void physics::gaugeObservables::writePlaqAndPolyToFile(int iter,  const std::string& filename)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile.width(8);
	outputToFile << iter << "\t";
	outputToFile.precision(15);
	outputToFile << plaquette << "\t" << plaquette_temporal << "\t" << plaquette_spatial << "\t" << polyakov.re << "\t" << polyakov.im << "\t" << sqrt(polyakov.re * polyakov.re + polyakov.im * polyakov.im) << std::endl;
	outputToFile.close();
}

void physics::gaugeObservables::measureTransportcoefficientKappa(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters)
{
	kappa = physics::algorithms::kappa_clover(gaugefield, parameters.get_beta());
	writeTransportcoefficientKappaToFile(parameters.get_transportcoefficientKappaFilename(), iteration);
}

void physics::gaugeObservables::writeTransportcoefficientKappaToFileUsingOpenOutputStream(int iteration)
{
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << iteration << "\t" << kappa << std::endl;
}

void physics::gaugeObservables::writeTransportcoefficientKappaToFile(std::string filename, int iteration)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if ( outputToFile.is_open() ) {
	  writeTransportcoefficientKappaToFileUsingOpenOutputStream(iteration);
	  outputToFile.close();
	} else {
		logger.warn() << "Could not open " << filename;
		File_Exception(filename.c_str());
	}
}

void physics::gaugeObservables::measureRectangles(physics::lattices::Gaugefield& gf, int iteration, meta::Inputparameters parameters)
{
	rectangles = gf.rectangles();
	writeRectanglesToFile(iteration, parameters.get_rectanglesFilename());
}

void physics::gaugeObservables::writeRectanglesToFile(int iter, const std::string& filename)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile.width(8);
	outputToFile << iter << "\t";
	outputToFile.precision(15);
	outputToFile << rectangles << std::endl;
	outputToFile.close();
}
