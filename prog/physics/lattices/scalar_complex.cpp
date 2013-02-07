/** @file
 * Implementation of some operations for physics::lattices::Scalar<hmc_complex>
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "scalar_complex.hpp"

#include <stdexcept>
#include "../../hardware/code/spinors.hpp"

//void physics::lattices::add(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right)
//{
//	auto dest_bufs = dest->get_buffers();
//	size_t num_bufs = dest_bufs.size();
//	auto left_bufs = left.get_buffers();
//	auto right_bufs = right.get_buffers();
//
//	if(num_bufs != left_bufs.size() || num_bufs != right_bufs.size()) {
//		throw std::invalid_argument("All arguments must use the same number of devices.");
//	}
//	for(size_t i = 0; i < num_bufs; ++i) {
//		auto code = dest_bufs[i]->get_device()->get_spinor_code();
//		code->product(dest_bufs[i], left_bufs[i], right_bufs[i]);
//	}
//}
//
//void physics::lattices::subtract(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& minuend, const Scalar<hmc_complex>& subtrahend);
//{
//	auto dest_bufs = dest->get_buffers();
//	size_t num_bufs = dest_bufs.size();
//	auto minuend_bufs = minuend.get_buffers();
//	auto subtrahend_bufs = subtrahend.get_buffers();
//
//	if(num_bufs != minuend_bufs.size() || num_bufs != subtrahend_bufs.size()) {
//		throw std::invalid_argument("All arguments must use the same number of devices.");
//	}
//	for(size_t i = 0; i < num_bufs; ++i) {
//		auto code = dest_bufs[i]->get_device()->get_spinor_code();
//		code->product(dest_bufs[i], minuend_bufs[i], subtrahend_bufs[i]);
//	}
//}

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
		auto code = dest_bufs[i]->get_device()->get_spinor_code();
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
		auto code = dest_bufs[i]->get_device()->get_spinor_code();
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
		auto code = dest_bufs[i]->get_device()->get_spinor_code();
		code->set_complex_to_float_device(src_bufs[i], dest_bufs[i]);
	}
}
