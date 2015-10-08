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

#ifndef _META_PARAMETERS_FERMION_HPP_
#define _META_PARAMETERS_FERMION_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersFermion {
public:
	common::action get_fermact() const noexcept;
	common::action get_fermact_mp() const noexcept;
	double get_kappa() const noexcept;
	double get_mass() const noexcept;
	double get_mu() const noexcept;
	double get_csw() const noexcept;
	double get_kappa_mp() const noexcept;
	double get_mu_mp() const noexcept;
	double get_csw_mp() const noexcept;
	double get_theta_fermion_spatial() const noexcept;
	double get_theta_fermion_temporal() const noexcept;
	double get_chem_pot_re() const noexcept;
	double get_chem_pot_im() const noexcept;
	bool get_use_chem_pot_re() const noexcept;
	bool get_use_chem_pot_im() const noexcept;
	bool get_use_eo() const noexcept;
	bool get_use_merge_kernels_fermion() const noexcept;
	bool get_use_merge_kernels_spinor() const noexcept;

private:
	po::options_description options;

	double kappa;
	double mass; //staggered quark mass
	double mu;
	double csw;
	double kappa_mp;
	double mu_mp;
	double csw_mp;
	double theta_fermion_spatial;
	double theta_fermion_temporal;
	bool use_chem_pot_re;
	bool use_chem_pot_im;
	double chem_pot_re;
	double chem_pot_im;
	bool use_eo;
	bool use_merge_kernels_fermion;
	bool use_merge_kernels_spinor;

protected:
	ParametersFermion();
	virtual ~ParametersFermion();
	ParametersFermion(ParametersFermion const&) = delete;
	ParametersFermion & operator=(ParametersFermion const&) = delete;
	po::options_description & getOptions();

	common::action fermact;
	common::action fermact_mp;
};

}

#endif
