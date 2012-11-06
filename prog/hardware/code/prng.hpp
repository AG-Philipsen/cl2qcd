/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_PRNG_
#define _HARDWARE_CODE_PRNG_

#include "opencl_module.hpp"

#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class PRNG : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~PRNG();

	/**
	 * Get cl_mem object rndarray
	 * @return rndarray
	 */
	const hardware::buffers::PRNGBuffer& get_prng_buffer() const noexcept;

	ClSourcePackage get_sources() const noexcept;

#ifdef USE_PRNG_RANLUX
	/**
	 * Initialize the state of the PRNG with the given seed.
	 */
	void initialize(const hardware::buffers::PRNGBuffer * buffer, cl_uint seed);
#endif /* USE_PRNG_RANLUX */

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
	PRNG(const meta::Inputparameters& params, hardware::Device * device);

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
#elif defined(USE_PRNG_RANLUX)
	cl_kernel init_kernel;
#endif // USE_PRNG_???
};

}

}

#endif // _HARDWARE_CODE_PRNG_
