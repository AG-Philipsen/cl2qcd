/** @file
 * Generic utility functions
 */

#include "util.hpp"

#include "../host_geometry.h"

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
size_t meta::get_spinorfieldsize(const Inputparameters& params)
{
	return get_vol4d(params);
}
size_t meta::get_eoprec_spinorfieldsize(const Inputparameters& params)
{
	return get_spinorfieldsize(params) / 2;
}

size_t meta::get_float_size(const Inputparameters& params)
{
	return params.get_precision() / 8;
}
size_t meta::get_mat_size(const Inputparameters& params)
{
	// TODO with rec12 this becomes 6
	return 9;
}

size_t meta::get_plaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * NDIM * (NDIM - 1) ) / 2.;
}
size_t meta::get_tplaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * (NDIM - 1));
}
size_t meta::get_splaq_norm(const Inputparameters& params)
{
	return (get_vol4d(params) * (NDIM - 1) * (NDIM - 2)) / 2. ;
}
size_t meta::get_rect_norm(const Inputparameters& params)
{
	return NDIM * (NDIM - 1) * get_vol4d(params);
}
size_t meta::get_poly_norm(const Inputparameters& params)
{
	return get_volspace(params);
}
size_t meta::get_flop_complex_mult() noexcept {
	return 6;
}
size_t meta::get_flop_su3_su3() noexcept {
	return (get_flop_complex_mult() * NC + (NC - 1) * 2) * NC * NC;
}
size_t meta::get_flop_su3trace() noexcept {
	return (NC - 1) * 2;
}
size_t meta::get_flop_su3_su3vec() noexcept {
	//  1 entry: NC * complex mults and NC-1 complex adds
	//  NC entries total
	return (get_flop_complex_mult() * NC + (NC - 1) * 2) * NC;
}
size_t meta::get_su3algebrasize() noexcept {
	return NC * NC - 1;
}
double meta::get_c0(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case Inputparameters::tlsym:
		case Inputparameters::iwasaki:
		case Inputparameters::dbw2:
			return 1. - 8. * get_c1(params);
		default:
			return 1.;
	}
}
double meta::get_c1(const Inputparameters& params)
{
	switch(params.get_gaugeact()) {
		case Inputparameters::tlsym:
			return -0.083333333;
		case Inputparameters::iwasaki:
			return -0.331;
		case Inputparameters::dbw2:
			return -1.4069;
		default:
			return 0.;
	}
}
double meta::get_xi_0(const Inputparameters& params)
{
	double aniso = params.get_xi();
	double beta = params.get_beta();
	double eta = (1.002503 * aniso * aniso * aniso + .39100 * aniso * aniso + 1.47130 * aniso - 0.19231) /
	             (aniso * aniso * aniso + 0.26287 * aniso * aniso + 1.59008 * aniso - 0.18224);
	return aniso / (1. + (1. - 1. / aniso) * eta / 6. * (1 - 0.55055 * 2 * NC / beta) / (1 - 0.77810 * 2 * NC / beta) * 2 * NC / beta );
}
size_t meta::get_flop_spinor_spinor() noexcept {
	//  NDIM * NC * complex_mult + ( NDIM * NC -1 ) complex adds
	return NDIM * NC * get_flop_complex_mult() + (NDIM * NC - 1) * 2;
}
size_t meta::get_flop_spinor_sqnorm() noexcept {
	//  NDIM * NC * 0.5 complex_mult + ( NDIM * NC -1 ) real adds
	return NDIM * NC * get_flop_complex_mult() * 0.5 + (NC * NDIM - 1);
}
