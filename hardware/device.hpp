/** @file
 * Declaration of the hardware::Device class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _HARDWARE_DEVICE_HPP_
#define _HARDWARE_DEVICE_HPP_

#include "device_info.hpp"
#include <map>
#include "hardwareParameters.hpp"
#include "opencl_compiler.hpp"
#include "profiling_data.hpp"
#include "../common_header_files/types.h"
#include "size_4.hpp"

class MemObjectAllocationTracer;

namespace hardware {

class SynchronizationEvent;
class Device;
class OpenClCode; //including the header of this class causes some trouble!

namespace buffers {
// forward declaration for friend relation
class Buffer;
class ProxyBufferCache;
hardware::SynchronizationEvent copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *dest_origin, const size_t *src_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, size_t src_row_pitch, size_t src_slice_pitch, const std::vector<hardware::SynchronizationEvent>& events);
}
namespace transfer {
class DGMAGhostBuffer;
}

namespace code {
// forward decleration to improve decoupling and speed up compilation
class Gaugefield;
class PRNG;
class Real;
class Complex;
class Spinors;
class Spinors_staggered;
class Fermions;
class Fermions_staggered;
class Correlator;
class Correlator_staggered;
class Heatbath;
class Kappa;
class Gaugemomentum;
class Molecular_Dynamics;
class Buffer;
}

class OptimizationError {
	// thrown to indicate that optimization failed
};

/**
 * Representation of an OpenCL device we are working on.
 *
 * Allows the querying and manipulation of the system and
 * its hardware.
 */
class Device : public DeviceInfo {

	friend hardware::buffers::Buffer;
	friend hardware::buffers::ProxyBufferCache;
	friend hardware::transfer::DGMAGhostBuffer;
	friend void print_profiling(Device *, const std::string&, int);
	friend cl_command_queue profiling_data_test_command_queue_helper(const Device * device);
	friend hardware::SynchronizationEvent hardware::buffers::copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *dest_origin, const size_t *src_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, size_t src_row_pitch, size_t src_slice_pitch, const std::vector<hardware::SynchronizationEvent>& events);
	friend MemObjectAllocationTracer;

public:
	/**
	 * Initialize an OpenCL device.
	 */
	Device(cl_context, cl_device_id, size_4 grid_pos, size_4 grid_size, const hardware::OpenClCode & builderIn, const hardware::HardwareParametersInterface & parametersIn);

	~Device();

	// non-copyable
	Device& operator=(const Device&) = delete;
	Device(const Device&) = delete;
	Device() = delete;

	/**
	 * Create a kernel from source files.
	 *
	 * Usage:
	 * @code
	 * cl_kernel dummy = create_kernel("dummy") << "dummy.cl";
	 * @endcode
	 */
	TmpClKernel create_kernel(const char * const kernelName, std::string buildOptions) const;

	/**
	 * Enqueue a kernel on the device using the default number of global threads
	 */
	void enqueue_kernel(const cl_kernel kernel);

	/**
	 * Enqueue a kernel on the device using the default number of local threads
	 */
	void enqueue_kernel(const cl_kernel kernel, size_t global_threads);

	/**
	 * Enqueue a kernel on the device using the default given threads specifications
	 */
	void enqueue_kernel(const cl_kernel kernel, size_t global_threads, size_t local_threads);

	void enqueue_marker(cl_event *) const;

	/**
	 * Enqueue a barrier, preventing new jobs starting on this device until the given event finished
	 */
	void enqueue_barrier(const hardware::SynchronizationEvent& event) const;
	void enqueue_barrier(const hardware::SynchronizationEvent& event1, const hardware::SynchronizationEvent& event2) const;

	/**
	 * Recommend a stride for a given number of elements of a non-basic type (called complete type).
	 *
	 * @todo rename fct. to recommendedStride
	 */
	size_t recommend_stride(size_t elementsOfCompleteType, size_t sizeOfBasicStorageElement, size_t numbeOfBasicStorageElementsInCompleteType) const;

	bool is_profiling_enabled() const noexcept;

	/**
	 * Make sure all commands have been sent to the device.
	 */
	void flush() const;

	/**
	 * Make sure all operations on this device are finished.
	 */
	void synchronize() const;

	/**
	 * Get the profiling data for all executions of the given kernel on this device.
	 */
	ProfilingData get_profiling_data(const cl_kernel& desiredKernel) noexcept;

	/**
	 * @todo: Rename all these getter fcts.
	 * e.g. get_gaugefield_code -> getGaugefieldCode
	 */
	/**
	 * Get access to the specific kernels on this device.
	 */
	const hardware::code::Gaugefield * get_gaugefield_code();
	const hardware::code::PRNG * get_prng_code();
	const hardware::code::Real * get_real_code();
	const hardware::code::Complex * get_complex_code();
	const hardware::code::Spinors * get_spinor_code();
	const hardware::code::Spinors_staggered * get_spinor_staggered_code();
	const hardware::code::Fermions * get_fermion_code();
	const hardware::code::Fermions_staggered * get_fermion_staggered_code();
	const hardware::code::Gaugemomentum * get_gaugemomentum_code();
	const hardware::code::Molecular_Dynamics * get_molecular_dynamics_code();
	const hardware::code::Correlator * get_correlator_code();
	const hardware::code::Correlator_staggered * get_correlator_staggered_code();
	const hardware::code::Heatbath * get_heatbath_code();
	const hardware::code::Kappa * get_kappa_code();
	/**
	 * TODO technicall this should only be used by stuff in the buffers package
	 */
	const hardware::code::Buffer * get_buffer_code();

	/**
	 *  TODO work over fct. names
	 * Get the position of the device inside the device grid.
	 */
	size_4 get_grid_pos() const;

	/**
	 * Get the size of the device grid.
	 */
	size_4 get_grid_size() const;
	size_4 get_local_lattice_size() const;
	unsigned get_halo_size() const;

	/**
	 * Get the size of the lattice in device memory.
	 */
	size_4 get_mem_lattice_size() const;

private:
	const cl_context context;
	const hardware::HardwareParametersInterface * hardwareParameters;
	const hardware::OpenClCode * openClCodeBuilder;
	cl_command_queue command_queue;

	/**
	 * Allow easy use of the command queue
	 */
	operator cl_command_queue() const noexcept;
	cl_command_queue get_queue() const noexcept;

	std::map<cl_kernel, ProfilingData> profiling_data;

	/**
	 * Pointers to specific code objects.
	 * Initialized on demand.
	 */
	const hardware::code::Gaugefield * gaugefield_code;
	const hardware::code::PRNG * prng_code;
	const hardware::code::Real * real_code;
	const hardware::code::Complex * complex_code;
	const hardware::code::Spinors * spinor_code;
	const hardware::code::Spinors_staggered * spinor_staggered_code;
	const hardware::code::Fermions * fermion_code;
	const hardware::code::Fermions_staggered * fermion_staggered_code;
	const hardware::code::Gaugemomentum * gaugemomentum_code;
	const hardware::code::Molecular_Dynamics * molecular_dynamics_code;
	const hardware::code::Correlator * correlator_code;
	const hardware::code::Correlator_staggered * correlator_staggered_code;
	const hardware::code::Heatbath * heatbath_code;
	const hardware::code::Kappa * kappa_code;
	const hardware::code::Buffer * buffer_code;

	/**
	 *  TODO work over member names
	 * The position of the device in the device grid.
	 */
	const size_4 grid_pos;

	/**
	 * The size of the device grid.
	 */
	const size_4 grid_size;
	const size_4 local_lattice_size;
	const unsigned halo_size;

	/**
	 * The size of the lattice in device memory.
	 */
	const size_4 mem_lattice_size;

	// memory usage tracing
	size_t allocated_bytes;
	size_t max_allocated_bytes;
	size_t allocated_hostptr_bytes;

	void markMemReleased(bool host, size_t size);
	void markMemAllocated(bool host, size_t size);
};

	/**
	 * Print the profiling information of kernels run on the given device.
	 */
	void print_profiling(Device * device, const std::string& filenameToWriteTo, int deviceId);
}

#endif /* _HARDWARE_DEVICE_HPP_ */
