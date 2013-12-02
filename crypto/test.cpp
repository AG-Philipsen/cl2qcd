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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Crypto
#include <boost/test/unit_test.hpp>

/**
 * Utility function to run a string through md5_buffer() and compare
 * the output with the given reference hash.
 *
 * \param string The string to calculate the MD5 of.
 * \param ref The reference MD5 value of the string.
 */
void testMD5(std::string string, std::string ref);

BOOST_AUTO_TEST_CASE( MD5 )
{
	/*
	 * MD5 Test Suite from RFC1321: http://ds.internic.net:/rfc/rfc1321.txt
	 *
	 * MD5 ("") = d41d8cd98f00b204e9800998ecf8427e
	 * MD5 ("a") = 0cc175b9c0f1b6a831c399e269772661
	 * MD5 ("abc") = 900150983cd24fb0d6963f7d28e17f72
	 * MD5 ("message digest") = f96b697d7cb7938d525a2f31aaf161d0
	 * MD5 ("abcdefghijklmnopqrstuvwxyz") = c3fcd3d76192e4007dfb496cca67e13b
	 * MD5 ("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789") =
	 * d174ab98d277d9f5a5611c2c9f419d9f
	 * MD5 ("123456789012345678901234567890123456789012345678901234567890123456
	 * 78901234567890") = 57edf4a22be3c955ac49da2e2107b67a
	 */
	testMD5("", "d41d8cd98f00b204e9800998ecf8427e");
	testMD5("a", "0cc175b9c0f1b6a831c399e269772661");
	testMD5("abc", "900150983cd24fb0d6963f7d28e17f72");
	testMD5("message digest", "f96b697d7cb7938d525a2f31aaf161d0");
	testMD5("abcdefghijklmnopqrstuvwxyz", "c3fcd3d76192e4007dfb496cca67e13b");
	testMD5("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789", "d174ab98d277d9f5a5611c2c9f419d9f");
	testMD5("12345678901234567890123456789012345678901234567890123456789012345678901234567890", "57edf4a22be3c955ac49da2e2107b67a");
}

void testMD5(std::string string, std::string ref)
{
	BOOST_CHECK_EQUAL(crypto::md5(string), ref);
}
