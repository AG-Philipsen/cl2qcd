/** @file
 * Implementation of explicit fermionamtrix operations
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "fermionmatrix.hpp"

void physics::fermionmatrix::M_wilson(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_wilson_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa);
	}
	out->update_halo();
}

void physics::fermionmatrix::M_tm_plus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_plus_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa, mubar);
	}

	out->update_halo();
}

void physics::fermionmatrix::M_tm_minus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_minus_device(in_bufs[i], out_bufs[i], gf_bufs[i], kappa, mubar);
	}

	out->update_halo();
}

void physics::fermionmatrix::M_tm_inverse_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_inverse_sitediagonal_device(in_bufs[i], out_bufs[i], mubar);
	}

	auto const valid_halo_width = in.get_valid_halo_width();
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::fermionmatrix::M_tm_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_sitediagonal_device(in_bufs[i], out_bufs[i], mubar);
	}

	auto const valid_halo_width = in.get_valid_halo_width();
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::fermionmatrix::M_tm_inverse_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_inverse_sitediagonal_minus_device(in_bufs[i], out_bufs[i], mubar);
	}

	auto const valid_halo_width = in.get_valid_halo_width();
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::fermionmatrix::M_tm_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar)
{
	auto out_bufs = out->get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->M_tm_sitediagonal_minus_device(in_bufs[i], out_bufs[i], mubar);
	}

	auto const valid_halo_width = in.get_valid_halo_width();
	if(valid_halo_width) {
		out->mark_halo_clean(valid_halo_width);
	} else {
		out->mark_halo_dirty();
	}
}

void physics::fermionmatrix::dslash(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, int evenodd, hmc_float kappa)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

#ifdef ASYNC_HALO_UPDATES
	auto update = in.require_halo_async(1);

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->dslash_eo_inner(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd, kappa);
	}

	update.finalize();

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->dslash_eo_boundary(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd, kappa);
	}
#else
	in.require_halo(1);

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->getFermionCode();
		fermion_code->dslash_eo_device(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd, kappa);
	}
#endif

	out->mark_halo_dirty();
}
