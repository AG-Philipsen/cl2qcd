/** @file
 * Definition of a counter
 */

#ifndef _META_COUNTER_
#define _META_COUNTER_

namespace meta {

	/**
	 * Generic counter
	 */
	class Counter {

	public:
		/**
		 * Create the counter, initialized to 0
		 */
		Counter() noexcept;

		/**
		 * Increment the counter
		 */
		Counter& operator+=(const unsigned&) noexcept;

		/**
		 * Increment the counter
		 */
		Counter& operator++() noexcept;

		/**
		 * Evaluation operator
		 */
		operator unsigned() const noexcept;

		/**
		 * Reset the counter
		 */
		void reset() noexcept;

	private:
		/**
		 * The actual value of the counter.
		 */
		unsigned value;
	};
}

#endif /* _META_COUNTER_ */
