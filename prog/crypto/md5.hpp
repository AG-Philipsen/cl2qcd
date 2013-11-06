#ifndef _CRYPTO_HPP_
#define _CRYPTO_HPP_

#include <string>

namespace crypto {

/**
 * Calculate the MD5 of the given string and return its hex representation.
 */
std::string md5(std::string);

}

#endif
