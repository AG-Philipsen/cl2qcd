#include "gaugeObservables.h"

void physics::gaugeObservables::writePlaqAndPolyToFile(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter,  const std::string& filename)
{
	std::ofstream outputFile(filename.c_str(), std::ios::app);
	if(!outputFile.is_open()) throw File_Exception(filename);
	outputFile.width(8);
	outputFile << iter << "\t";
	outputFile.precision(15);
	outputFile << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	outputFile.close();
}

void physics::gaugeObservables::measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params)
{
	hmc_float plaq;
	hmc_float tplaq;
	hmc_float splaq;
	hmc_complex pol;

	gf.gaugeobservables(&plaq, &tplaq, &splaq, &pol);

	if ( params.get_print_to_screen() ){
	  logger.info() << iter << '\t' << plaq << '\t' << tplaq << '\t' << splaq << '\t' << pol.re << '\t' << pol.im << '\t' << sqrt(pol.re * pol.re + pol.im * pol.im);
	}
	writePlaqAndPolyToFile(plaq, tplaq, splaq, pol, iter, filename);
}

void physics::gaugeObservables::measurePlaqAndPoly(physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params)
{
	const std::string filename = meta::get_gauge_obs_file_name(params,  "");
	measurePlaqAndPoly(gf, iter, filename, params);
}

void physics::gaugeObservables::measureGaugeObservables(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters)
{
  measurePlaqAndPoly(gaugefield, iteration, parameters);
  if ( parameters.get_measure_transportcoefficient_kappa() ) {
    measureTransportcoefficientKappa(gaugefield, iteration, parameters);
  }
}


void physics::gaugeObservables::measureTransportcoefficientKappa(physics::lattices::Gaugefield& gaugefield, int iteration, meta::Inputparameters parameters)
{
	double kappa = 0;
	kappa = physics::algorithms::kappa_clover(gaugefield, parameters.get_beta());
	writeTransportcoefficientKappaToFile(kappa, "kappa_clover.dat", iteration);
}

void physics::gaugeObservables::writeTransportcoefficientKappaToFileUsingOpenOutputStream(hmc_float kappa, int iteration)
{
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << iteration << "\t" << kappa << std::endl;
}

void physics::gaugeObservables::writeTransportcoefficientKappaToFile(hmc_float kappa, std::string filename, int iteration)
{
	outputToFile.open(filename.c_str(), std::ios::app);
	if ( outputToFile.is_open() ) {
	  writeTransportcoefficientKappaToFileUsingOpenOutputStream(kappa, iteration);
		outputToFile.close();
	} else {
		logger.warn() << "Could not open " << filename;
		File_Exception(filename.c_str());
	}
}

