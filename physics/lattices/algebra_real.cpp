/** @file
 * Implementation of some algebra operations for physics::lattices::Scalar<hmc_float>
 * together with physics::lattices::Vector<hmc_float>
 *
 * Copyright (c) 2014 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "algebra_real.hpp"

#include <stdexcept>
#include "../../hardware/device.hpp"
#include "../../hardware/code/real.hpp"


void physics::lattices::update_zeta_cgm(const Vector<hmc_float>* out, const Vector<hmc_float>& zeta_prev, const Vector<hmc_float>& zeta_prev_prev, const Scalar<hmc_float>& sbeta_prev, const Scalar<hmc_float>& sbeta_pres, const Scalar<hmc_float>& salpha_prev, const Vector<hmc_float>& sigma, const int numeq)
{
	auto out_bufs = out->get_buffers();
	size_t num_bufs = out_bufs.size();
	auto zeta_prev_bufs = zeta_prev.get_buffers();
	auto zeta_prev_prev_bufs = zeta_prev_prev.get_buffers();
	auto sbeta_prev_bufs = sbeta_prev.get_buffers();
	auto sbeta_pres_bufs = sbeta_pres.get_buffers();
	auto salpha_prev_bufs = salpha_prev.get_buffers();
	auto sigma_bufs = sigma.get_buffers();

	if(num_bufs != zeta_prev_bufs.size() || num_bufs != zeta_prev_prev_bufs.size() ||
	   num_bufs != sbeta_prev_bufs.size() || num_bufs != sbeta_pres_bufs.size() ||
	   num_bufs != salpha_prev_bufs.size() || num_bufs != sigma_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = out_bufs[i]->get_device()->get_real_code();
		code->update_zeta_cgm_device(zeta_prev_bufs[i], zeta_prev_prev_bufs[i], sbeta_prev_bufs[i], sbeta_pres_bufs[i], salpha_prev_bufs[i], sigma_bufs[i], numeq, out_bufs[i]);
	}
}


void physics::lattices::update_beta_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& sbeta_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& zeta_prev, const int numeq)
{
	auto out_bufs = out->get_buffers();
	size_t num_bufs = out_bufs.size();
	auto sbeta_pres_bufs = sbeta_pres.get_buffers();
	auto zeta_pres_bufs = zeta_pres.get_buffers();
	auto zeta_prev_bufs = zeta_prev.get_buffers();

	if(num_bufs != zeta_pres_bufs.size() || num_bufs != zeta_prev_bufs.size() ||
	   num_bufs != sbeta_pres_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = out_bufs[i]->get_device()->get_real_code();
		code->update_beta_cgm_device(sbeta_pres_bufs[i], zeta_pres_bufs[i], zeta_prev_bufs[i], numeq, out_bufs[i]);
	}
}


void physics::lattices::update_alpha_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& salpha_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& beta_pres, const Vector<hmc_float>& zeta_prev, const Scalar<hmc_float>& sbeta_pres, const int numeq)
{
	auto out_bufs = out->get_buffers();
	size_t num_bufs = out_bufs.size();
	auto salpha_pres_bufs = salpha_pres.get_buffers();
	auto zeta_pres_bufs = zeta_pres.get_buffers();
	auto beta_pres_bufs = beta_pres.get_buffers();
	auto zeta_prev_bufs = zeta_prev.get_buffers();
	auto sbeta_pres_bufs = sbeta_pres.get_buffers();
	
	if(num_bufs != salpha_pres_bufs.size() || num_bufs != zeta_pres_bufs.size() ||
	   num_bufs != beta_pres_bufs.size() || num_bufs != zeta_prev_bufs.size() ||
	   num_bufs != sbeta_pres_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = out_bufs[i]->get_device()->get_real_code();
		code->update_alpha_cgm_device(salpha_pres_bufs[i], zeta_pres_bufs[i], beta_pres_bufs[i], zeta_prev_bufs[i], sbeta_pres_bufs[i], numeq, out_bufs[i]);
	}
}

