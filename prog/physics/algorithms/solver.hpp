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

#include "../Fermionmatrix.hpp"

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
	Solver(Fermionmatrix & M, cl_mem x, cl_mem b, cl_mem gf, hmc_float prec, int max_iter, hardware::Device device);
	
	virtual ~Solver();
	
	/**
	 * Solve the system
	 */
	virtual int solve();
	
	/**
	 * Access to private members
	 */
	size_t get_bytes() const noexcept;
	cl_mem get_x() const noexcept;
	cl_mem get_b() const noexcept;
	cl_mem get_gf() const noexcept;
	Matrix_Function * get_M() const noexcept;
	hmc_float get_prec() const noexcept;
	int get_iter();
	int get_iter_max() const noexcept;
	
	/**
	 * Get the device this buffer is located on
	 */
	hardware::Device * get_device() const noexcept;
	
      private:
	/**
	 * The ingredients to the linear system to be solved
	 */
	const cl_mem x;
	const cl_mem b;
	const cl_mem gf;
	const Fermionmatrix & M;
	
	/**
	 * The acceptance precision for the solve
	 */
	const hmc_float prec;
	
	/**
	 * Counter for the iterations during the solver
	 */
	int iter;
	
	/**
	 * Maximal number of iterations
	 */
	const int max_iter;
	
	/**
	 * The OpenCL device the solver is working on.
	 */
	hardware::Device * const device;
	
      protected:
      };
    }
  }
}  
#endif /* _PHYSICS_ALGORITHMS_ */
