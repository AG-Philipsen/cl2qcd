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

void physics::lattices::access_real_vector_element(const Scalar<hmc_float>* out, const Vector<hmc_float>& in, const int index)
{
	auto out_bufs = out->get_buffers();
	size_t num_bufs = out_bufs.size();
	auto in_bufs = in.get_buffers();

	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = out_bufs[i]->get_device()->getRealCode();
		code->set_real_to_vector_element_device(in_bufs[i], index, out_bufs[i]);
	}
}  

void physics::lattices::access_real_vector_element(const Vector<hmc_float>* out, const Scalar<hmc_float>& in, const int index)
{
	auto out_bufs = out->get_buffers();
	size_t num_bufs = out_bufs.size();
	auto in_bufs = in.get_buffers();

	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = out_bufs[i]->get_device()->getRealCode();
		code->set_vector_element_to_real_device(in_bufs[i], index, out_bufs[i]);
	}
}  

void physics::lattices::add(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& left, const Scalar<hmc_float>& right)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto left_bufs = left.get_buffers();
	auto right_bufs = right.get_buffers();

	if(num_bufs != left_bufs.size() || num_bufs != right_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getRealCode();
		code->set_real_to_sum_device(left_bufs[i], right_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::subtract(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& minuend, const Scalar<hmc_float>& subtrahend)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto minuend_bufs = minuend.get_buffers();
	auto subtrahend_bufs = subtrahend.get_buffers();

	if(num_bufs != minuend_bufs.size() || num_bufs != subtrahend_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getRealCode();
		code->set_real_to_difference_device(minuend_bufs[i], subtrahend_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::multiply(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& left, const Scalar<hmc_float>& right)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto left_bufs = left.get_buffers();
	auto right_bufs = right.get_buffers();

	if(num_bufs != left_bufs.size() || num_bufs != right_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getRealCode();
		code->set_real_to_product_device(left_bufs[i], right_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::divide(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& numerator, const Scalar<hmc_float>& denominator)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto numerator_bufs = numerator.get_buffers();
	auto denominator_bufs = denominator.get_buffers();

	if(num_bufs != numerator_bufs.size() || num_bufs != denominator_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getRealCode();
		code->set_real_to_ratio_device(numerator_bufs[i], denominator_bufs[i], dest_bufs[i]);
	}
}


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
		auto code = out_bufs[i]->get_device()->getRealCode();
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
		auto code = out_bufs[i]->get_device()->getRealCode();
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
		auto code = out_bufs[i]->get_device()->getRealCode();
		code->update_alpha_cgm_device(salpha_pres_bufs[i], zeta_pres_bufs[i], beta_pres_bufs[i], zeta_prev_bufs[i], sbeta_pres_bufs[i], numeq, out_bufs[i]);
	}
}

size_t physics::lattices::get_flops_update_cgm(const std::string quantity, const int Neqs, const hardware::System& system)
{
	// assert single system
	auto devices = system.get_devices();
	auto real_code = devices[0]->getRealCode();
	if(quantity == "alpha")
	  return real_code->get_flop_size_update("update_alpha_cgm", Neqs);
	else if(quantity == "beta")
	  return real_code->get_flop_size_update("update_beta_cgm", Neqs);
	else if(quantity == "zeta")
	  return real_code->get_flop_size_update("update_zeta_cgm", Neqs);
	else
	  throw Invalid_Parameters("Quantity unknown in get_flops_update_cgm", "alpha, beta or zeta", quantity);
}



