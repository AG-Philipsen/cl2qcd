/** @file
 * Implementation of the transfer method factory
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "transfer.hpp"
#include "ocl_copy.hpp"
#include "async_ocl_copy.hpp"
#include "dgma.hpp"
#include "../../logger.hpp"

std::unique_ptr<hardware::Transfer> hardware::create_transfer(hardware::Device * const src, hardware::Device * const dest, hardware::System const & system)
{
// required headers for DirectGMA might be missing
// in that case it is unavailable
#ifdef CL_MEM_BUS_ADDRESSABLE_AMD
	try {
		return std::unique_ptr<hardware::Transfer>(new transfer::DirectGMA(src, dest, system));
	} catch(hardware::transfer::DGMAUnsupported) {
		logger.warn() << "DirectGMA is unavailable. Falling back to standard transfer methods.";
	}
#endif
#ifdef ASYNC_HALO_UPDATES
	return std::unique_ptr<hardware::Transfer>(new transfer::AsyncOclCopy(src, dest, system));
#else
	return std::unique_ptr<hardware::Transfer>(new transfer::OclCopy(src, dest));
#endif
}
