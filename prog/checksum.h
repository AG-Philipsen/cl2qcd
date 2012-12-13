/** @file
 * Dekleration of the checksum class
 */

#ifndef _CHECKSUM_H_
#define _CHECKSUM_H_

#include "types.h"

#include <boost/crc.hpp>
#include <ostream>

class Checksum {

public:
	/**
	 * Create class given precalculated checksum values
	 */
	Checksum(uint32_t suma, uint32_t sumb) noexcept : suma(suma), sumb(sumb) { };

	/**
	 * Create an empty checksum (checksum of 0 bytes)
	 */
	Checksum() noexcept : suma(0), sumb(0) { };

	inline uint32_t get_suma() const noexcept {
		return suma;
	}

	inline uint32_t get_sumb() const noexcept {
		return sumb;
	}

	/**
	 * Add another site to the mix
	 *
	 * @param rank
	 *        tmlqcd uses -> ((t*LZ+z)*LY+y)*LX+x
	 */
	inline void accumulate(const char* buf, size_t size, uint32_t rank) {
		uint32_t rank29 = rank % 29;
		uint32_t rank31 = rank % 31;

		boost::crc_32_type worker;
		worker.process_bytes(buf, size);
		uint32_t work = worker.checksum();

		suma ^= work<<rank29 | work>>(32-rank29);
		sumb ^= work<<rank31 | work>>(32-rank31);
	}

	bool operator==(const Checksum& other) {
		return suma == other.suma && sumb == other.sumb;
	}

	bool operator!=(const Checksum& other) {
		return suma != other.suma || sumb != other.sumb;
	}

private:
	uint32_t suma;
	uint32_t sumb;
};

inline std::ostream& operator<<(std::ostream& out, const Checksum& sum) noexcept {
	out << "SumA: " << std::hex << sum.get_suma() << ", SumB: " << std::hex << sum.get_sumb();
	return out;
}

#endif /* _CHECKSUM_H_ */
