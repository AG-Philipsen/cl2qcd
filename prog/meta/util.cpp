/** @file
 * Generic utility functions
 */

#include "util.hpp"

size_t meta::get_volspace(const Inputparameters& params)
{
	size_t ns = params.get_nspace();
	return ns * ns * ns;
}
size_t meta::get_vol4d(const Inputparameters& params)
{
	return get_volspace(params) * params.get_ntime();
}
bool meta::get_use_rectangles(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case meta::Inputparameters::action::tlsym:
		case meta::Inputparameters::action::iwasaki:
		case meta::Inputparameters::action::dbw2:
			return true;
		default:
			return false;
	}
}
hmc_float meta::get_mubar(const Inputparameters& params)
{
	return 2. * params.get_kappa() * params.get_mu();
}
hmc_float meta::get_mubar_mp(const Inputparameters& params)
{
	return 2. * params.get_kappa_mp() * params.get_mu_mp();
}
