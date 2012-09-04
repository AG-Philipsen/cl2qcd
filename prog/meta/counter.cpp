#include "counter.hpp"

meta::Counter::Counter() noexcept
:
value(0)
{
	// already initialized
}

meta::Counter& meta::Counter::operator+=(const unsigned& inc) noexcept {
	value += inc;
	return *this;
}

meta::Counter& meta::Counter::operator++() noexcept {
	value++;
	return *this;
}

meta::Counter::operator unsigned() const noexcept
{
	return value;
}

void meta::Counter::reset() noexcept {
	value = 0;
}
