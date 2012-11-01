#ifndef _FERMIONMATRIXH_
#define _FERMIONMATRIXH_

#include "../opencl_module.h"

#include "../hardware/buffers/plain.hpp"
#include "../hardware/buffers/su3.hpp"
#include "../hardware/buffers/spinor.hpp"
#include "../host_use_timer.h"


/**
 * this is the definition of the class "Fermionmatrix"
 */
class Fermionmatrix {
protected:
	Opencl_Module * that;

	Fermionmatrix(Opencl_Module * that) : that(that) { };

public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;

	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_Flops() const = 0;

	/**
	 * Get the net bytes read / written by this function.
	 */
	virtual cl_ulong get_Bytes() const = 0;
};

#endif /* _FERMIONMATRIXH_ */
