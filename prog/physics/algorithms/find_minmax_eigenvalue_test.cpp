/** @file
 * Implementation of the algorithm to find min and max eigenvalues of an operator
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::find_minmax_eigenvalue
#include <boost/test/unit_test.hpp>
