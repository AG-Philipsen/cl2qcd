/*
 * Copyright 2013 Matthias Bach
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

#include "md5.hpp"

#include <iomanip>
#include <sstream>
#include <type_traits>
#include <nettle/md5.h>

std::string crypto::md5(std::string plain)
{
	struct md5_ctx ctx;
	md5_init(&ctx);
	static_assert(std::is_same<unsigned char, uint8_t>::value, "char should be the same as uint8_t!");
	md5_update(&ctx, plain.length(), reinterpret_cast<const uint8_t *>(plain.c_str()));
	uint8_t res[MD5_DIGEST_SIZE];
	md5_digest(&ctx, MD5_DIGEST_SIZE, res);

	std::ostringstream fingerprint;
	for(uint8_t b: res) {
		fingerprint << std::setfill('0') << std::setw(2) << std::hex << static_cast<uint32_t>(b);
	}
	return fingerprint.str();
}
