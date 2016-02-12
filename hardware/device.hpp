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
#include "../geometry/latticeGrid.hpp"

class MemObjectAllocationTracer;

namespace hardware {

class SynchronizationEvent;
class Device;
class OpenClCode; //including the header of this class causes some trouble!

namespace buffers {
// forward declaration for friend relation
class Buffer;
class ProxyBufferCache;
hardware::SynchronizationEvent copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *destOrigin, const size_t *srcOrigin, const size_t *region, size_t destRowPitch, size_t destSlicePitch, size_t srcRowPitch, size_t srcSlicePitch, const std::vector<hardware::SynchronizationEvent>& events);
}
namespace transfer {
class DGMAGhostBuffer;
}

namespace code {
// forward decleration to improve decoupling and speed up compilation
class Gaugefield;
class Prng;
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
	friend void printProfiling(Device *, const std::string&, int);
	friend cl_command_queue profilingDataTestCommandQueueHelper(const Device * device);
	friend hardware::SynchronizationEvent hardware::buffers::copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *destOrigin, const size_t *srcOrigin, const size_t *region, size_t destRowPitch, size_t destSlicePitch, size_t srcRowPitch, size_t srcSlicePitch, const std::vector<hardware::SynchronizationEvent>& events);
	friend MemObjectAllocationTracer;

public:
	/**
	 * Initialize an OpenCL device.
	 */
	Device(cl_context, cl_device_id, LatticeGridIndex lI, LatticeGrid lG, const hardware::OpenClCode & builderIn, const hardware::HardwareParametersInterface & parametersIn);

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
	TmpClKernel createKernel(const char * const kernelName, std::string buildOptions) const;

	/**
	 * Enqueue a kernel on the device using the default number of global threads
	 */
	void enqueueKernel(const cl_kernel kernel) const;

	/**
	 * Enqueue a kernel on the device using the default number of local threads
	 */
	void enqueueKernel(const cl_kernel kernel, size_t globalThreads) const;

	/**
	 * Enqueue a kernel on the device using the default given threads specifications
	 */
	void enqueue_kernel(const cl_kernel kernel, size_t globalThreads, size_t localThreads) const;

	void enqueueMarker(cl_event *) const;

	/**
	 * Enqueue a barrier, preventing new jobs starting on this device until the given event finished
	 */
	void enqueueBarrier(const hardware::SynchronizationEvent& event) const;
	void enqueue_barrier(const hardware::SynchronizationEvent& event1, const hardware::SynchronizationEvent& event2) const;

	/**
	 * Recommend a stride for a given number of elements of a non-basic type (called complete type).
	 *
	 * @todo rename fct. to recommendedStride
	 */
	size_t recommendStride(size_t elementsOfCompleteType, size_t sizeOfBasicStorageElement, size_t numbeOfBasicStorageElementsInCompleteType) const;

	bool isProfilingEnabled() const noexcept;

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
	ProfilingData getProfilingData(const cl_kernel& desiredKernel) const noexcept;

	/**
	 * @todo: Rename all these getter fcts.
	 * e.g. get_gaugefield_code -> getGaugefieldCode
	 */
	/**
	 * Get access to the specific kernels on this device.
	 */
	const hardware::code::Gaugefield * getGaugefieldCode() const;
	const hardware::code::Prng * getPrngCode() const;
	const hardware::code::Real * getRealCode() const;
	const hardware::code::Complex * getComplexCode() const;
	const hardware::code::Spinors * getSpinorCode() const;
	const hardware::code::Spinors_staggered * getSpinorStaggeredCode() const;
	const hardware::code::Fermions * getFermionCode() const;
	const hardware::code::Fermions_staggered * getFermionStaggeredCode() const;
	const hardware::code::Gaugemomentum * getGaugemomentumCode() const;
	const hardware::code::Molecular_Dynamics * getMolecularDynamicsCode() const;
	const hardware::code::Correlator * getCorrelatorCode() const;
	const hardware::code::Correlator_staggered * getCorrelatorStaggeredCode() const;
	const hardware::code::Heatbath * getHeatbathCode() const;
	const hardware::code::Kappa * getKappaCode() const;
	/**
	 * TODO technically this should only be used by stuff in the buffers package
	 */
	const hardware::code::Buffer * getBufferCode() const;

	/**
	 *  TODO work over fct. names
	 * Get the position of the device inside the device grid.
	 */
	LatticeGridIndex getGridPos() const;

	/**
	 * Get the size of the device grid.
	 * @todo: move to system!
	 */
	LatticeGrid getGridSize() const;
	LocalLatticeExtents getLocalLatticeExtents() const;
	unsigned getHaloExtent() const;
	LocalLatticeMemoryExtents getLocalLatticeMemoryExtents() const;

private:
	const hardware::OpenClCode * openClCodeBuilder;
	const hardware::HardwareParametersInterface * hardwareParameters;
	const cl_context context;
	cl_command_queue command_queue;

	/**
	 * Allow easy use of the command queue
	 */
	operator cl_command_queue() const noexcept;
	cl_command_queue get_queue() const noexcept;

	mutable std::map<cl_kernel, ProfilingData> profiling_data;

	/**
	 * Pointers to specific code objects.
	 * Initialized on demand.
	 */
	mutable hardware::code::Gaugefield * gaugefield_code;
	mutable hardware::code::Prng * prng_code;
	mutable hardware::code::Real * real_code;
	mutable hardware::code::Complex * complex_code;
	mutable hardware::code::Spinors * spinor_code;
	mutable hardware::code::Spinors_staggered * spinor_staggered_code;
	mutable hardware::code::Fermions * fermion_code;
	mutable hardware::code::Fermions_staggered * fermion_staggered_code;
	mutable hardware::code::Gaugemomentum * gaugemomentum_code;
	mutable hardware::code::Molecular_Dynamics * molecular_dynamics_code;
	mutable hardware::code::Correlator * correlator_code;
	mutable hardware::code::Correlator_staggered * correlator_staggered_code;
	mutable hardware::code::Heatbath * heatbath_code;
	mutable hardware::code::Kappa * kappa_code;
	mutable hardware::code::Buffer * buffer_code;

	/**
	 *  TODO work over member names
	 * The position of the device in the device grid.
	 */
	LatticeGridIndex latticeGridIndex;

	/**
	 * The size of the device grid.
	 */
	const LatticeGrid latticeGridExtents;
	const LocalLatticeExtents localLatticeExtents;
	const unsigned haloSize;

	/**
	 * The size of the lattice in device memory.
	 */
	const LocalLatticeMemoryExtents localLatticeMemoryExtents;

	// memory usage tracing
	mutable size_t allocated_bytes;
	mutable size_t max_allocated_bytes;
	mutable size_t allocated_hostptr_bytes;

	void markMemReleased(bool host, size_t size) const;
	void markMemAllocated(bool host, size_t size) const;
};

	/**
	 * Print the profiling information of kernels run on the given device.
	 */
	void printProfiling(Device * device, const std::string& filenameToWriteTo, int deviceId);
}

#endif /* _HARDWARE_DEVICE_HPP_ */
