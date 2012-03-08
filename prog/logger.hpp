#ifndef _LOGGERHPP_
#define _LOGGERHPP_
/**
 * @file
 *
 * Singleton for logging to screen
 */

#include "einhard/einhard.hpp"

/**
 * Singleton logger instance for logging to stdout.
 */
extern einhard::Logger<einhard::LOG_LEVEL> logger;


#endif //_LOGGERHPP_
