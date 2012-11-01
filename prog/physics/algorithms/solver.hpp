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
#include "../../logger.hpp"

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
	Solver(hmc_float prec, int _max_iter, int _iter_refresh) : acc_prec(prec), iter_max(_max_iter), iter_refresh(_iter_refresh), iter(0) {};
	
	virtual ~Solver();
	
	/**
	 * Access to private members
	 */
	hmc_float get_prec() const noexcept;
	int get_iter();
	int get_iter_max() const noexcept;
	int get_iter_refresh() const noexcept;

	/**
	 * Solve the system
	 */
	virtual bool solve() const;
	
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
	const int iter_max;

	/**
	 * Number of iterations where the solution guess is refreshed
	 */
	const int iter_refresh;
	
      protected:
      };
      /**
       * A Solver to solve system without even-odd preconditioning
       */
      class Solver_noneo : public Solver{
      public:
	Solver_noneo(const physics::fermionmatrix::Fermionmatrix & _f, const hardware::buffers::Plain<spinor> * _x, const hardware::buffers::Plain<spinor> * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver(prec, max_iter, iter_refresh), x(_x), b(_b), gf(_gf), f(_f) {};
	/**
	 * Access to private members
	 */
	hardware::buffers::Plain<spinor> * get_x() const noexcept;
	hardware::buffers::Plain<spinor> * get_b() const noexcept;
	hardware::buffers::SU3 * get_gf() const noexcept;
	physics::fermionmatrix::Fermionmatrix * get_M() const noexcept;
      protected:
	/**
	 * The ingredients to the linear system to be solved
	 */
	const hardware::buffers::Plain<spinor> * x;
	const hardware::buffers::Plain<spinor> * b;
	const hardware::buffers::SU3 * gf;
	const physics::fermionmatrix::Fermionmatrix & f;
      };
      /**
       * A Solver to solve system using even-odd preconditioning
       */
      class Solver_eo : public Solver{
      public:
	Solver_eo(const physics::fermionmatrix::Fermionmatrix_eo & _f, const hardware::buffers::Spinor * _x, const hardware::buffers::Spinor * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver(prec, max_iter, iter_refresh), x(_x), b(_b), gf(_gf), f(_f) {};
	/**
	 * Access to private members
	 */
	hardware::buffers::Spinor * get_x() const noexcept;
	hardware::buffers::Spinor * get_b() const noexcept;
	hardware::buffers::SU3 * get_gf() const noexcept;
	physics::fermionmatrix::Fermionmatrix * get_M() const noexcept;
      protected:
	/**
	 * The ingredients to the linear system to be solved
	 */
	const hardware::buffers::Spinor * x;
	const hardware::buffers::Spinor * b;
	const hardware::buffers::SU3 * gf;
	const physics::fermionmatrix::Fermionmatrix_eo & f;
      };
      /**
       * The Conjugate Gradient (CG) solver
       */
      class Cg : public Solver_noneo{
	Cg(const physics::fermionmatrix::Fermionmatrix & _f, const hardware::buffers::Plain<spinor> * _x, const hardware::buffers::Plain<spinor> * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver_noneo(_f, _x, _b, _gf, prec, max_iter, iter_refresh) {};
      };
      class Cg_eo : public Solver_eo{
      public:
	Cg_eo(const physics::fermionmatrix::Fermionmatrix_eo & _f, const hardware::buffers::Spinor * _x, const hardware::buffers::Spinor * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver_eo(_f, _x, _b, _gf, prec, max_iter, iter_refresh) {};
	bool solve() const override;
      };
      /**
       * The Biconjugate Gradient Stabilized (BiCGStab) solver
       */
      class Bicgstab : public Solver_noneo{
	Bicgstab(const physics::fermionmatrix::Fermionmatrix & _f, const hardware::buffers::Plain<spinor> * _x, const hardware::buffers::Plain<spinor> * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver_noneo(_f, _x, _b, _gf, prec, max_iter, iter_refresh) {};
      public:
      };
      class Bicgstab_eo : public Solver_eo{
	Bicgstab_eo(const physics::fermionmatrix::Fermionmatrix_eo & _f, const hardware::buffers::Spinor * _x, const hardware::buffers::Spinor * _b, const hardware::buffers::SU3 * _gf, hmc_float prec, int max_iter, int iter_refresh)
	  : Solver_eo(_f, _x, _b, _gf, prec, max_iter, iter_refresh) {};
      };
    }
  }
}  
#endif /* _PHYSICS_ALGORITHMS_ */
