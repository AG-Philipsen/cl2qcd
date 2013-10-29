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
