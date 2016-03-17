/*
 * Copyright 2016 Francesca Cuteri
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

#include "gaugefield.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../host_functionality/host_operations_gaugefield.h"
#include "../../host_functionality/host_random.h"
#include "../device.hpp"
#include "../buffers/halo_update.hpp"
#include "../code/gaugefield.hpp"
#include "../../geometry/parallelization.hpp"

static Matrixsu3 random_matrixsu3();

hardware::lattices::Gaugefield::Gaugefield(const hardware::System& system):
system(system), buffers(allocate_buffers()),unsmeared_buffers()
{}

hardware::lattices::Gaugefield::~Gaugefield()
{
	release_buffers(&buffers);
	release_buffers(&unsmeared_buffers);
}

const std::vector<const hardware::buffers::SU3 *> hardware::lattices::Gaugefield::get_buffers() const noexcept
{
	return buffers;
}

std::vector<const hardware::buffers::SU3 *> hardware::lattices::Gaugefield::allocate_buffers()
{
	using hardware::buffers::SU3;

	std::vector<const SU3 *> buffers;

	auto const devices = system.get_devices();
	for(auto device: devices) 
	{
		buffers.push_back(new SU3(device->getLocalLatticeMemoryExtents().getLatticeVolume() * 4, device)); //todo: do not calculate here!
	}
	return buffers;
}

void hardware::lattices::Gaugefield::release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers)
{
	for(auto buffer: *buffers) 
	{
		delete buffer;
	}
	buffers->clear();
}

void hardware::lattices::Gaugefield::send_gaugefield_to_buffers(const Matrixsu3 * const gf_host)
{
	logger.trace() << "importing gaugefield";
// 	if(buffers.size() == 1) {
// 		auto device = buffers[0]->get_device();
// 		device->getGaugefieldCode()->importGaugefield(buffers[0], gf_host);
// 		device->synchronize();
// 	} else {
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();
			TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrixsu3), device->getHaloExtent());

			if(buffers.size() == 1) device->getGaugefieldCode()->importGaugefield(buffer, gf_host);
			else{
				Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];
//				//todo: put these calls into own fct.! With smart pointers?
				memcpy(&mem_host[tmp2.getMainPartIndex_destination()]  , &gf_host[tmp2.getMainPartIndex_source()]  , tmp2.getMainPartSizeInBytes());
				memcpy(&mem_host[tmp2.getFirstHaloIndex_destination()] , &gf_host[tmp2.getFirstHaloPartIndex_source()] , tmp2.getHaloPartSizeInBytes());
				memcpy(&mem_host[tmp2.getSecondHaloIndex_destination()], &gf_host[tmp2.getSecondHaloPartIndex_source()], tmp2.getHaloPartSizeInBytes());

				device->getGaugefieldCode()->importGaugefield(buffer, mem_host);
				delete[] mem_host;
			}
			device->synchronize();
		}
// 	}
	logger.trace() << "import complete";
}

void hardware::lattices::Gaugefield::fetch_gaugefield_from_buffers(Matrixsu3 * const gf_host)
{
	logger.trace() << "fetching gaugefield";
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->getGaugefieldCode()->exportGaugefield(gf_host, buffers[0]);
		device->synchronize();
	} else {
		for(auto const buffer: buffers) {
			// fetch local part for each device
			auto device = buffer->get_device();
			Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];

			device->getGaugefieldCode()->exportGaugefield(mem_host, buffer);

			TemporalParallelizationHandlerLink tmp2(device->getGridPos(), device->getLocalLatticeExtents(), sizeof(Matrixsu3), device->getHaloExtent());
			memcpy(&gf_host[tmp2.getMainPartIndex_source()]  , &mem_host[tmp2.getMainPartIndex_destination()]  , tmp2.getMainPartSizeInBytes());

			delete[] mem_host;
		}
	}
}

void hardware::lattices::Gaugefield::update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<Matrixsu3>(buffers, system, NDIM);
}

void hardware::lattices::Gaugefield::update_halo_soa(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<Matrixsu3>(buffers, system, .5, 2 * NDIM);
}

void hardware::lattices::Gaugefield::update_halo() const
{
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(buffers, system);
		} else {
			update_halo_aos(buffers, system);
		}
	}
}

void hardware::lattices::Gaugefield::set_cold() const
{
	using hardware::Device;

	for(auto buffer: buffers) 
	{
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_cold(tmp, elems);
		const Device * device = buffer->get_device();
		device->getGaugefieldCode()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

void hardware::lattices::Gaugefield::set_cold(Matrixsu3 * field, size_t elems) const
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void hardware::lattices::Gaugefield::set_hot() const
{
	using hardware::Device;

	for(auto buffer: buffers) 
	{
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_hot(tmp, elems);
		const Device * device = buffer->get_device();
		device->getGaugefieldCode()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

void hardware::lattices::Gaugefield::set_hot(Matrixsu3 * field, size_t elems) const
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = random_matrixsu3();
	}
}

static Matrixsu3 random_matrixsu3()
{
	//this is like project_su3 from operations_gaugefield.cl, only that you do it with gaussian su3 vectors
	//taken also from the tmlqcd code
	Matrixsu3 out;
	hmc_complex a[NC];
	hmc_complex b[NC];
	hmc_complex c[NC];

	//get gaussian vector
	gaussianComplexVector(a, NC, 1.);
	//get norm 1
	hmc_float norm = 0;
	for(int i = 0; i < NC; i++) {
		norm += a[i].re * a[i].re + a[i].im * a[i].im;
	}
	for(int i = 0; i < NC; i++) {
		a[i].re /= sqrt(norm);
		a[i].im /= sqrt(norm);
	}

	for(;;) {
		//get another gaussian vector
		gaussianComplexVector(b, NC, 1.);
		//get norm 1
		norm = 0;
		for(int i = 0; i < NC; i++) {
			norm += b[i].re * b[i].re + b[i].im * b[i].im;
		}
		for(int i = 0; i < NC; i++) {
			b[i].re /= sqrt(norm);
			b[i].im /= sqrt(norm);
		}

		//calculate a* times b
		hmc_complex tmp = {0., 0.};
		for(int i = 0; i < NC; i++) {
			tmp.re += a[i].re * b[i].re + a[i].im * b[i].im;
			tmp.im += - a[i].im * b[i].re + a[i].re * b[i].im;
		}

		//project: b = tmp * a
		b[0].re -= (tmp.re * a[0].re - tmp.im * a[0].im);
		b[0].im -= (tmp.re * a[0].im + tmp.im * a[0].re);
		b[1].re -= (tmp.re * a[1].re - tmp.im * a[1].im);
		b[1].im -= (tmp.re * a[1].im + tmp.im * a[1].re);
		b[2].re -= (tmp.re * a[2].re - tmp.im * a[2].im);
		b[2].im -= (tmp.re * a[2].im + tmp.im * a[2].re);

		//get norm
		norm = 0;
		for(int i = 0; i < NC; i++) {
			norm += b[i].re * b[i].re + b[i].im * b[i].im;
		}
		norm = sqrt(norm);

		//check if norm is zero
		if (1.0 != (1.0 + norm))
			break;
	}

	// b = b/norm
	hmc_float fact = 1.0 / norm;
	b[0].re = fact * b[0].re;
	b[0].im = fact * b[0].im;
	b[1].re = fact * b[1].re;
	b[1].im = fact * b[1].im;
	b[2].re = fact * b[2].re;
	b[2].im = fact * b[2].im;

	//3rd vector, this is reconstruct12!!
	c[0].re = (a[1].re * b[2].re - a[1].im * b[2].im)
	          - (a[2].re * b[1].re - a[2].im * b[1].im);
	c[0].im = -(a[1].re * b[2].im + a[1].im * b[2].re)
	          + (a[2].re * b[1].im + a[2].im * b[1].re);

	c[1].re = (a[2].re * b[0].re - a[2].im * b[0].im)
	          - (a[0].re * b[2].re - a[0].im * b[2].im);
	c[1].im = -(a[2].re * b[0].im + a[2].im * b[0].re)
	          + (a[0].re * b[2].im + a[0].im * b[2].re);

	c[2].re = (a[0].re * b[1].re - a[0].im * b[1].im)
	          - (a[1].re * b[0].re - a[1].im * b[0].im);
	c[2].im = -(a[0].re * b[1].im + a[0].im * b[1].re)
	          + (a[1].re * b[0].im + a[1].im * b[0].re);

	//set values of matrix
	out.e00.re = a[0].re;
	out.e00.im = a[0].im;
	out.e01.re = a[1].re;
	out.e01.im = a[1].im;
	out.e02.re = a[2].re;
	out.e02.im = a[2].im;

	out.e10.re = b[0].re;
	out.e10.im = b[0].im;
	out.e11.re = b[1].re;
	out.e11.im = b[1].im;
	out.e12.re = b[2].re;
	out.e12.im = b[2].im;

	out.e20.re = c[0].re;
	out.e20.im = c[0].im;
	out.e21.re = c[1].re;
	out.e21.im = c[1].im;
	out.e22.re = c[2].re;
	out.e22.im = c[2].im;

	return out;
}

void hardware::lattices::Gaugefield::smear(unsigned int smearingSteps)
{
	unsmeared_buffers = allocate_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		auto device = buf->get_device();
		auto gf_code = device->getGaugefieldCode(); // should be like: get_gaugefield_code( HardwareParameters_gaugefield )

		hardware::buffers::copyData(unsmeared_buffers[i], buf);

		int rho_iter = smearingSteps;
		logger.debug() << "\t\tperform " << rho_iter << " steps of stout-smearing to the gaugefield...";

		//one needs a temporary gf to apply the smearing to
		const hardware::buffers::SU3 gf_tmp(buf->get_elements(), device);
		for(int i = 0; i < rho_iter - 1; i += 2) {
			gf_code->stout_smear_device(buf, &gf_tmp);
			gf_code->stout_smear_device(&gf_tmp, buf);
		}
		//if rho_iter is odd one has to copy ones more
		if(rho_iter % 2 == 1) {
			gf_code->stout_smear_device(buf, &gf_tmp);
			hardware::buffers::copyData(buf, &gf_tmp);
		}
	}
}

void hardware::lattices::Gaugefield::unsmear()
{
	if(unsmeared_buffers.size() == 0) {
		logger.warn() << "Tried to unsmear gaugefield that is not smeared.";
		return;
	}

	unsmeared_buffers = allocate_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		hardware::buffers::copyData(buf, unsmeared_buffers[i]);
	}

	release_buffers(&unsmeared_buffers);
}