/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULERANH_
#define _OPENCLMODULERANH_

#include "opencl_module.h"

#include "hardware/buffers/prng_buffer.hpp"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Ran : public Opencl_Module {
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
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string&) const { return 0; };

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string&) const { return 0; };

private:

	/**
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Ran(const meta::Inputparameters& params, hardware::Device * device);

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
