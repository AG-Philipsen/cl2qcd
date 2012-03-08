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

/**
 * Set the log level to the one described by the string.
 *
 * @arg log_level A string describing the log level. Not case sensetive, should be one of
 *                all, trace, debug, info, warn, error, fatal, off
 */
void switchLogLevel(const std::string log_level);

#endif //_LOGGERHPP_
