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
#include "../meta/inputparameters.hpp"
#include "hardwareParameters.hpp"
#include "opencl_compiler.hpp"
#include "profiling_data.hpp"
#include "../common_header_files/types.h"
#include "../meta/size_4.hpp"

class MemObjectAllocationTracer;

namespace hardware {

class SynchronizationEvent;
class Device;
class OpenClCode;

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
	 * Initialize the device.
	 *
	 * \param context context for the device
	 * \param device_id id of the device to initialize
	 * \param inputparams the input parameters of the application
	 * \param enable_profiling enable profiling on this device
	 */
	Device(cl_context, cl_device_id, size_4 grid_pos, size_4 grid_size, const meta::Inputparameters&, bool enable_profiling = false);

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
	 *
	 * @param kernel_name The name of the kernel to create.
	 * @param build_opts Build options for the kernel
	 */
	TmpClKernel create_kernel(const char * const kernel_name, std::string build_opts) const;

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

	/**
	 * Enqueue a kernel on the device using the default number of global threads
	 */
	void enqueue_marker(cl_event *) const;

	/**
	 * Enqueue a barrier, preventing new jobs starting on this device until the given event finished
	 */
	void enqueue_barrier(const hardware::SynchronizationEvent& event) const;
	void enqueue_barrier(const hardware::SynchronizationEvent& event1, const hardware::SynchronizationEvent& event2) const;

	/**
	 * Recommend a stride for the given number of elements of the given type
	 *
	 * \param elems      The number of elements to be stored
	 * \param type_size  The size of the basic storage element
	 * \param lane_count How many of the basic storage elements make up the complete type
	 * \return The recommended stride in elements
	 */
	size_t recommend_stride(size_t elems, size_t type_size, size_t lane_count) const;

	/**
	 * Query whether profiling is enabled on this device.
	 */
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
	 *
	 * \param kernel The kernel for which to get the profiling data
	 * \return Profiling data for the given kernel on this device
	 */
	ProfilingData get_profiling_data(const cl_kernel& kernel) noexcept;

	/**
	 * Get access to the gaugefield kernels on this device.
	 */
	const hardware::code::Gaugefield * get_gaugefield_code();

	/**
	 * Get access to the prng kernels on this device.
	 */
	const hardware::code::PRNG * get_prng_code();

	/**
	 * Get access to the real kernels on this device.
	 */
	const hardware::code::Real * get_real_code();
	
	/**
	 * Get access to the complex kernels on this device.
	 */
	const hardware::code::Complex * get_complex_code();

	/**
	 * Get access to the spinor kernels on this device.
	 */
	const hardware::code::Spinors * get_spinor_code();

	/**
	 * Get access to the staggered spinor kernels on this device.
	 */
	const hardware::code::Spinors_staggered * get_spinor_staggered_code();

	/**
	 * Get access to the fermion kernels on this device.
	 */
	const hardware::code::Fermions * get_fermion_code();

	/**
	 * Get access to the staggered fermion kernels on this device.
	 */
	const hardware::code::Fermions_staggered * get_fermion_staggered_code();

	/**
	 * Get access to the gaugemomentum kernels on this device.
	 */
	const hardware::code::Gaugemomentum * get_gaugemomentum_code();

	/**
	 * Get access to the molecular dynamics kernels on this device.
	 */
	const hardware::code::Molecular_Dynamics * get_molecular_dynamics_code();

	/**
	 * Get access to the correlator kernels on this device.
	 */
	const hardware::code::Correlator * get_correlator_code();
	
	/**
	 * Get access to the correlator_staggered kernels on this device.
	 */
	const hardware::code::Correlator_staggered * get_correlator_staggered_code();

	/**
	 * Get access to the heatbath kernels on this device.
	 */
	const hardware::code::Heatbath * get_heatbath_code();

	/**
	 * Get access to the kappa kernels on this device.
	 */
	const hardware::code::Kappa * get_kappa_code();

	/**
	 * Get access to the buffer kernels on this device.
	 *
	 * TODO technicall this should only be used by stuff in the buffers package
	 */
	const hardware::code::Buffer * get_buffer_code();

	/**
	 * Get the position of the device inside the device grid.
	 */
	size_4 get_grid_pos() const;

	/**
	 * Get the size of the device grid.
	 */
	size_4 get_grid_size() const;

	/**
	 * Get the size of the local lattice.
	 */
	size_4 get_local_lattice_size() const;

	/**
	 * Get the size of the halo
	 */
	unsigned get_halo_size() const;

	/**
	 * Get the size of the lattice in device memory.
	 */
	size_4 get_mem_lattice_size() const;

private:
	/**
	 * The OpenCL context to be used by this device.
	 */
	const cl_context context;

	/**
	 * The input parameters of the application.
	 */
	const meta::Inputparameters& params;
	const hardware::HardwareParametersInterface * hardwareParameters;
	const hardware::OpenClCode * openClCodeBuilder;

	/**
	 * The command queue used to perform operations on this device
	 */
	cl_command_queue command_queue;

	/**
	 * Whether profiling is enabled
	 */
	const bool profiling_enabled;

	/**
	 * Allow easy use of the command queue
	 */
	operator cl_command_queue() const noexcept;

	/**
	 * Allow easy use of the command queue
	 */
	cl_command_queue get_queue() const noexcept;

	/**
	 * Kernel profiling data
	 */
	std::map<cl_kernel, ProfilingData> profiling_data;

	/**
	 * Pointer to the gaugefield code.
	 * Initialized on demand.
	 */
	const hardware::code::Gaugefield * gaugefield_code;

	/**
	 * Pointer to the prng code.
	 * Initialized on demand.
	 */
	const hardware::code::PRNG * prng_code;

	/**
	 * Pointer to the real code.
	 * Initialized on demand.
	 */
	const hardware::code::Real * real_code;
	
	/**
	 * Pointer to the complex code.
	 * Initialized on demand.
	 */
	const hardware::code::Complex * complex_code;

	/**
	 * Pointer to the spinor code.
	 * Initialized on demand.
	 */
	const hardware::code::Spinors * spinor_code;

	/**
	 * Pointer to the staggered spinor code.
	 * Initialized on demand.
	 */
	const hardware::code::Spinors_staggered * spinor_staggered_code;

	/**
	 * Pointer to the fermion code.
	 * Initialized on demand.
	 */
	const hardware::code::Fermions * fermion_code;

	/**
	 * Pointer to the staggered fermion code.
	 * Initialized on demand.
	 */
	const hardware::code::Fermions_staggered * fermion_staggered_code;

	/**
	 * Pointer to the gaugemomentum code.
	 * Initialized on demand.
	 */
	const hardware::code::Gaugemomentum * gaugemomentum_code;

  	/**
	 * Pointer to the molecular dynamics code.
	 * Initialized on demand.
	 */
	const hardware::code::Molecular_Dynamics * molecular_dynamics_code;

	/**
	 * Pointer to the correlator code.
	 * Initialized on demand.
	 */
	const hardware::code::Correlator * correlator_code;
	
	/**
	 * Pointer to the staggered correlator code.
	 * Initialized on demand.
	 */
	const hardware::code::Correlator_staggered * correlator_staggered_code;

	/**
	 * Pointer to the heatbath code.
	 * Initialized on demand.
	 */
	const hardware::code::Heatbath * heatbath_code;

	/**
	 * Pointer to the kappa code.
	 * Initialized on demand.
	 */
	const hardware::code::Kappa * kappa_code;

	/**
	 * Pointer to the buffer code
	 * Initialized on demand.
	 */
	const hardware::code::Buffer * buffer_code;

	/**
	 * The position of the device in the device grid.
	 */
	const size_4 grid_pos;

	/**
	 * The size of the device grid.
	 */
	const size_4 grid_size;

	/**
	 * The size of the local lattice.
	 */
	const size_4 local_lattice_size;

	/**
	 * Get the size of the halo
	 */
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
	 *
	 * \param device The device the kernels ran on
	 * \param filename The file to write the profiling information to
	 * \param id The id to identify this device by
	 */
	void print_profiling(Device * device, const std::string& filename, int id);
}

#endif /* _HARDWARE_DEVICE_HPP_ */
