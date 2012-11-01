/** @file
 * Declaration of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_
#define _PHYSICS_ALGORITHMS_

#include "../../hardware/device.hpp"

#include <stdexcept>
#include "../../hardware/system.hpp"

#include "../fermionmatrix.hpp"

namespace physics {
  namespace algorithms {
    /**
     * this namespace contains methods to solve linear systems like
     * A * x = b
     * for x
     */
    namespace solver{
      
      /**
       * A generic Solver.
       */
      class Solver {
	
      public:
	Solver(hmc_float prec, int _max_iter) : acc_prec(prec), max_iter(_max_iter), iter(0) {};
	
	virtual ~Solver();
	
	/**
	 * Solve the system
	 */
	virtual int solve();
	
	/**
	 * Access to private members
	 */
	hmc_float get_prec() const noexcept;
	int get_iter();
	int get_iter_max() const noexcept;
	
      private:
	/**
	 * The acceptance precision for the solve
	 */
	const hmc_float acc_prec;
	
	/**
	 * Counter for the iterations during the solver
	 */
	int iter;
	
	/**
	 * Maximal number of iterations
	 */
	const int max_iter;
	
      protected:
      };
      /**
       * A Solver to solve system without even-odd preconditioning
       */
      class Solver_noneo : public Solver{
      public:
	Solver_noneo(physics::fermionmatrix::Fermionmatrix & _f, cl_mem _x, cl_mem _b, cl_mem _gf, hmc_float prec, int max_iter) : Solver(prec, max_iter), x(_x), b(_b), gf(_gf), f(_f) {};
	/**
	 * Access to private members
	 */
	cl_mem get_x() const noexcept;
	cl_mem get_b() const noexcept;
	cl_mem get_gf() const noexcept;
	physics::fermionmatrix::Fermionmatrix * get_M() const noexcept;
      private:
	/**
	 * The ingredients to the linear system to be solved
	 */
	const cl_mem x;
	const cl_mem b;
	const cl_mem gf;
	const physics::fermionmatrix::Fermionmatrix & f;
	
	
      };
    }
  }
}  
#endif /* _PHYSICS_ALGORITHMS_ */
