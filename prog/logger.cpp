/**
 * @file
 *
 * Utility file to generate singleton logging object.
 */

#include "logger.hpp"

#include <string>
#include <algorithm>
#include <cctype>
#include <map>

einhard::Logger<einhard::LOG_LEVEL> logger(einhard::ALL);

std::string toLower(std::string in);

void switchLogLevel(const std::string log_level)
{
	std::map<std::string, einhard::LogLevel> log_level_map;
	for(einhard::LogLevel level = einhard::ALL; level <= einhard::OFF; level = static_cast<einhard::LogLevel>(level + 1)) {
		log_level_map[toLower(einhard::getLogLevelString(level))] = level;
	}
	logger.setVerbosity(log_level_map[toLower(log_level)]);
}

std::string toLower(std::string in)
{
	std::transform(in.begin(), in.end(), in.begin(), ::tolower);
	return in;
}
