/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULERANH_
#define _OPENCLMODULERANH_

#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "opencl_compiler.hpp"

#include "opencl_module_gaugefield.h"

#include "exceptions.h"

#include "hardware/buffers/prng_buffer.hpp"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Ran : public Opencl_Module_Gaugefield {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Ran();

	/**
	 * Get cl_mem object rndarray
	 * @return rndarray
	 */
	const hardware::buffers::PRNGBuffer& get_prng_buffer() const noexcept;

	ClSourcePackage get_sources() const noexcept;

protected:
	/**
	 * @param[in] params points to an instance of inputparameters
	 *
	 * @fixme Needs to be private
	 */
	Opencl_Module_Ran(const meta::Inputparameters& params, hardware::Device * device);
private:

	/**
	 * A set of sources required to use the PRNG.
	 */
	ClSourcePackage prng_code;

	const hardware::buffers::PRNGBuffer prng_buffer;

#ifdef USE_PRNG_NR3
	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndstate_to_device(nr3_state_dev* host_rndarray) const;

	/**
	 * Copy the RNG state from the OpenCL buffer.
	 *
	 * @param[out] rndarray The RNG copy target
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndstate_from_device(nr3_state_dev* rndarray) const;

	nr3_state_dev* rndarray;
	size_t sizeof_rndarray;
#endif // USE_PRNG_NR3
};

#endif //OPENCLMODULERANH
