/** @file
 * Generic utility functions
 */

#ifndef _META_UTIL_
#define _META_UTIL_

#include "inputparameters.hpp"

namespace meta {
	size_t get_volspace(const Inputparameters&);
	size_t get_vol4d(const Inputparameters&);
	bool get_use_rectangles(const Inputparameters& params);
}

#endif /* META_UTIL_ */
