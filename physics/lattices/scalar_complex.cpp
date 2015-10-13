/** @file
 * Implementation of some operations for physics::lattices::Scalar<hmc_complex>
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "scalar_complex.hpp"

#include <stdexcept>
#include "../../hardware/device.hpp"
#include "../../hardware/code/complex.hpp"

void physics::lattices::add(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto left_bufs = left.get_buffers();
	auto right_bufs = right.get_buffers();

	if(num_bufs != left_bufs.size() || num_bufs != right_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getComplexCode();
		code->set_complex_to_sum_device(left_bufs[i], right_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::subtract(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& minuend, const Scalar<hmc_complex>& subtrahend)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto minuend_bufs = minuend.get_buffers();
	auto subtrahend_bufs = subtrahend.get_buffers();

	if(num_bufs != minuend_bufs.size() || num_bufs != subtrahend_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getComplexCode();
		code->set_complex_to_difference_device(minuend_bufs[i], subtrahend_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::multiply(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto left_bufs = left.get_buffers();
	auto right_bufs = right.get_buffers();

	if(num_bufs != left_bufs.size() || num_bufs != right_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getComplexCode();
		code->set_complex_to_product_device(left_bufs[i], right_bufs[i], dest_bufs[i]);
	}
}
void physics::lattices::divide(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& numerator, const Scalar<hmc_complex>& denominator)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto numerator_bufs = numerator.get_buffers();
	auto denominator_bufs = denominator.get_buffers();

	if(num_bufs != numerator_bufs.size() || num_bufs != denominator_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getComplexCode();
		code->set_complex_to_ratio_device(numerator_bufs[i], denominator_bufs[i], dest_bufs[i]);
	}
}

void physics::lattices::convert(const Scalar<hmc_complex>* dest, const Scalar<hmc_float>& src)
{
	auto dest_bufs = dest->get_buffers();
	size_t num_bufs = dest_bufs.size();
	auto src_bufs = src.get_buffers();

	if(num_bufs != src_bufs.size()) {
		throw std::invalid_argument("All arguments must use the same number of devices.");
	}
	for(size_t i = 0; i < num_bufs; ++i) {
		auto code = dest_bufs[i]->get_device()->getComplexCode();
		code->set_complex_to_float_device(src_bufs[i], dest_bufs[i]);
	}
}
