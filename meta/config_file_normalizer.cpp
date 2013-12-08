/** @file
 * Implementation of utilities to normalize input file contents
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <boost/regex.hpp>

#include "config_file_normalizer.hpp"

/**
 * Load the contents of the given inputstream into a string.
 */
static std::string load_stream(std::istream& is);

/**
 * Build the regular expression and replacement string for the given map of aliases
 */
static std::pair<boost::regex, std::string> build_regex(std::map<std::string, std::string>);

void meta::ConfigFileNormalizer::add_alias(std::string alias, std::string mapped)
{
	aliases[alias] = mapped;
}

std::string meta::ConfigFileNormalizer::operator() (std::string file) const
{
	std::ifstream config_file(file.c_str());
	if(!config_file.is_open()) {
		throw std::invalid_argument(file + " could not be opened.");
	}
	return this->operator()(config_file);
}

std::string meta::ConfigFileNormalizer::operator() (std::istream& input) const
{
	std::string orig_config_file = load_stream(input);

	std::ostringstream normalized_file;
	std::ostream_iterator<char, char> normalized_file_sink(normalized_file);

	auto regex = build_regex(aliases);
	boost::regex_replace(normalized_file_sink, orig_config_file.begin(), orig_config_file.end(), regex.first, regex.second, boost::match_default | boost::format_all);

	return normalized_file.str();
}

static std::string load_stream(std::istream& is)
{
	std::string s;
	s.reserve(is.rdbuf()->in_avail());
	char c;
	while(is.get(c)) {
		if(s.capacity() == s.size())
			s.reserve(s.capacity() * 3);
		s.append(1, c);
	}
	return s;
}

static std::pair<boost::regex, std::string> build_regex(std::map<std::string, std::string> aliases)
{
	std::ostringstream regex;
	std::ostringstream replacements;

	int index = 0;
for(auto alias: aliases) {
		++index; // index goes from 1 to n
		if(index > 1) {
			regex << '|';
		}

		regex << "(^\\s*" << alias.first << "\\s*=)";
		replacements << "(?" << index << alias.second << "=)";
	}
	return std::make_pair<>(boost::regex(regex.str()), replacements.str());
}
