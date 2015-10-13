/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include "parametersFermion.hpp"

bool meta::ParametersFermion::get_use_chem_pot_re() const noexcept
{
	return use_chem_pot_re;
}
bool meta::ParametersFermion::get_use_chem_pot_im() const noexcept
{
	return use_chem_pot_im;
}

common::action meta::ParametersFermion::get_fermact() const noexcept
{
	return fermact;
}
common::action meta::ParametersFermion::get_fermact_mp() const noexcept
{
	return fermact_mp;
}
double meta::ParametersFermion::get_kappa() const noexcept
{
	return kappa;
}
double meta::ParametersFermion::get_mass() const noexcept
{
	return mass;
}
double meta::ParametersFermion::get_mu() const noexcept
{
	return mu;
}
double meta::ParametersFermion::get_csw() const noexcept
{
	return csw;
}
double meta::ParametersFermion::get_kappa_mp() const noexcept
{
	return kappa_mp;
}
double meta::ParametersFermion::get_mu_mp() const noexcept
{
	return mu_mp;
}
double meta::ParametersFermion::get_csw_mp() const noexcept
{
	return csw_mp;
}

double meta::ParametersFermion::get_theta_fermion_spatial() const noexcept
{
	return theta_fermion_spatial;
}
double meta::ParametersFermion::get_theta_fermion_temporal() const noexcept
{
	return theta_fermion_temporal;
}
double meta::ParametersFermion::get_chem_pot_re() const noexcept
{
	return chem_pot_re;
}
double meta::ParametersFermion::get_chem_pot_im() const noexcept
{
	return chem_pot_im;
}
bool meta::ParametersFermion::get_use_eo() const noexcept
{
	return use_eo;
}

bool meta::ParametersFermion::get_use_merge_kernels_fermion() const noexcept
{
	return use_merge_kernels_fermion;
}
bool meta::ParametersFermion::get_use_merge_kernels_spinor() const noexcept
{
	return use_merge_kernels_spinor;
}

meta::ParametersFermion::ParametersFermion()
	: options("Fermion options")
{
	options.add_options()
	("fermact", po::value<std::string>()->default_value("wilson"))
	("fermact_mp", po::value<std::string>()->default_value("wilson"))
	//todo: change this default value!
	("kappa", po::value<double>(&kappa)->default_value(0.125))
	("mass", po::value<double>(&mass)->default_value(0.1))
	("mu", po::value<double>(&mu)->default_value(0.006))
	("csw", po::value<double>(&csw)->default_value(0.))
	("kappa_mp", po::value<double>(&kappa_mp)->default_value(0.125))
	("mu_mp", po::value<double>(&mu_mp)->default_value(0.006))
	("csw_mp", po::value<double>(&csw_mp)->default_value(0.))
	("theta_fermion_spatial", po::value<double>(&theta_fermion_spatial)->default_value(0.))
	("theta_fermion_temporal", po::value<double>(&theta_fermion_temporal)->default_value(0.))
	("chem_pot_re", po::value<double>(&chem_pot_re)->default_value(0.))
	("chem_pot_im", po::value<double>(&chem_pot_im)->default_value(0.))
	("use_chem_pot_re", po::value<bool>(&use_chem_pot_re)->default_value(false))
	("use_chem_pot_im", po::value<bool>(&use_chem_pot_im)->default_value(false))
	("use_eo", po::value<bool>(&use_eo)->default_value(true))
	("use_merge_kernels_spinor", po::value<bool>(&use_merge_kernels_spinor)->default_value(false), "Use kernel merging for spinor kernels")
	("use_merge_kernels_fermion", po::value<bool>(&use_merge_kernels_fermion)->default_value(false), "Use kernel merging for fermion kernels");
}

meta::ParametersFermion::~ParametersFermion() = default;

po::options_description & meta::ParametersFermion::getOptions()
{
	return options;
}
