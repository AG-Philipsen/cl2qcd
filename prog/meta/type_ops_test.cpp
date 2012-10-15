/** @file
 * Testcases for the type operations
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "type_ops.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE type_ops
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(hmc_complex_ops)
{
	const hmc_complex one = {1., 0.};
	const hmc_complex zero = {0., 0.};
	const hmc_complex ei = {0., 1.};

	BOOST_REQUIRE_EQUAL(one, one);
	BOOST_REQUIRE_NE(one, zero);
	BOOST_REQUIRE_NE(one, ei);
	BOOST_REQUIRE_NE(zero, one);
	BOOST_REQUIRE_EQUAL(zero, zero);
	BOOST_REQUIRE_NE(zero, ei);
	BOOST_REQUIRE_NE(ei, one);
	BOOST_REQUIRE_NE(ei, zero);
	BOOST_REQUIRE_EQUAL(ei, ei);
}

BOOST_AUTO_TEST_CASE(Matrixsu3_ops)
{
	const hmc_complex c_one = {1., 0.};
	const hmc_complex c_zero = {0., 0.};
	const hmc_complex c_ei = {0., 1.};

	const Matrixsu3 one = {c_one, c_zero, c_zero, c_zero, c_one, c_zero, c_zero, c_zero, c_one};
	const Matrixsu3 zero = {c_zero, c_zero, c_zero, c_zero, c_zero, c_zero, c_zero, c_zero, c_zero};
	const Matrixsu3 ei = {c_one, c_ei, c_zero, c_zero, c_one, c_zero, c_zero, c_zero, c_one};

	BOOST_REQUIRE_EQUAL(one, one);
	BOOST_REQUIRE_NE(one, zero);
	BOOST_REQUIRE_NE(one, ei);
	BOOST_REQUIRE_NE(zero, one);
	BOOST_REQUIRE_EQUAL(zero, zero);
	BOOST_REQUIRE_NE(zero, ei);
	BOOST_REQUIRE_NE(ei, one);
	BOOST_REQUIRE_NE(ei, zero);
	BOOST_REQUIRE_EQUAL(ei, ei);
}

BOOST_AUTO_TEST_CASE(su3vec_ops)
{
	const hmc_complex c_one = {1., 0.};
	const hmc_complex c_zero = {0., 0.};
	const hmc_complex c_ei = {0., 1.};

	const su3vec one = {c_one, c_zero, c_zero};
	const su3vec zero = {c_zero, c_zero, c_zero};
	const su3vec ei = {c_one, c_ei, c_zero};

	BOOST_REQUIRE_EQUAL(one, one);
	BOOST_REQUIRE_NE(one, zero);
	BOOST_REQUIRE_NE(one, ei);
	BOOST_REQUIRE_NE(zero, one);
	BOOST_REQUIRE_EQUAL(zero, zero);
	BOOST_REQUIRE_NE(zero, ei);
	BOOST_REQUIRE_NE(ei, one);
	BOOST_REQUIRE_NE(ei, zero);
	BOOST_REQUIRE_EQUAL(ei, ei);
}

BOOST_AUTO_TEST_CASE(spinor_ops)
{
	const hmc_complex c_one = {1., 0.};
	const hmc_complex c_zero = {0., 0.};
	const hmc_complex c_ei = {0., 1.};

	const su3vec v_one = {c_one, c_zero, c_zero};
	const su3vec v_zero = {c_zero, c_zero, c_zero};
	const su3vec v_ei = {c_one, c_ei, c_zero};

	const spinor one = {v_one, v_one, v_one, v_one};
	const spinor zero = {v_zero, v_zero, v_zero, v_zero};
	const spinor ei = {v_one, v_ei, v_zero, v_ei};

	BOOST_REQUIRE_EQUAL(one, one);
	BOOST_REQUIRE_NE(one, zero);
	BOOST_REQUIRE_NE(one, ei);
	BOOST_REQUIRE_NE(zero, one);
	BOOST_REQUIRE_EQUAL(zero, zero);
	BOOST_REQUIRE_NE(zero, ei);
	BOOST_REQUIRE_NE(ei, one);
	BOOST_REQUIRE_NE(ei, zero);
	BOOST_REQUIRE_EQUAL(ei, ei);
}

BOOST_AUTO_TEST_CASE(ae_ops)
{
	const hmc_float c_one = {1.};
	const hmc_float c_zero = {0.};
	const hmc_float c_ei = { -1.};

	const ae one = {c_one, c_zero, c_zero, c_zero, c_one, c_zero, c_zero, c_zero};
	const ae zero = {c_zero, c_zero, c_zero, c_zero, c_zero, c_zero, c_zero, c_zero};
	const ae ei = {c_one, c_ei, c_zero, c_zero, c_one, c_zero, c_zero, c_zero};

	BOOST_REQUIRE_EQUAL(one, one);
	BOOST_REQUIRE_NE(one, zero);
	BOOST_REQUIRE_NE(one, ei);
	BOOST_REQUIRE_NE(zero, one);
	BOOST_REQUIRE_EQUAL(zero, zero);
	BOOST_REQUIRE_NE(zero, ei);
	BOOST_REQUIRE_NE(ei, one);
	BOOST_REQUIRE_NE(ei, zero);
	BOOST_REQUIRE_EQUAL(ei, ei);
}
