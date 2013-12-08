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

void physics::gaugeObservables::measurePlaqAndPoly(const physics::lattices::Gaugefield& gf, int iter, const std::string& filename, meta::Inputparameters params)
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

void physics::gaugeObservables::measurePlaqAndPoly(const physics::lattices::Gaugefield& gf, int iter, meta::Inputparameters params)
{
	const std::string filename = meta::get_gauge_obs_file_name(params,  "");
	measurePlaqAndPoly(gf, iter, filename, params);
}
