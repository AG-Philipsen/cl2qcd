/**
 * @file
 *
 * Utility file to generate singleton logging object.
 */

#include "logger.hpp"

einhard::Logger<einhard::LOG_LEVEL> logger(einhard::ALL);
