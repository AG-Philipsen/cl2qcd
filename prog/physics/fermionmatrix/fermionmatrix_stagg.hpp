/** @file
 * Definition of staggered fermionmatrix classes.
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#ifndef _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_STAGG_
#define _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_STAGG_

#include "../../hardware/code/fermions.hpp"

#include "../../host_use_timer.h"
#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/gaugefield.hpp"
#include "../../hardware/device.hpp"

/**
 * This is the definition of the classes "Fermionmatrix_stagg"
 */
namespace physics {

namespace fermionmatrix {

/*
 * Explicit Fermion operations
 */






/**
 * A generic staggered fermion matrix
 */
class Fermionmatrix_stagg_basic {
public:
	/**
	 * Get if the matrix is hermitian
	 */
	bool is_hermitian() const noexcept;
	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_flops() const = 0;

protected:
	Fermionmatrix_stagg_basic(const hardware::System& system, bool herm, hmc_float _mass = ARG_DEF) : _is_hermitian(herm), mass(_mass), system(system) { };

	/**
	 * Get the mass of the fermion.
	 */
	hmc_float get_mass() const noexcept;
	
	/**
	 * Get the system to operate on.
	 */
	const hardware::System& get_system() const noexcept;

private:
	/**
	 * Shows if matrix is hermitian.
	 */
	const bool _is_hermitian;

	/**
	 * Mass of the fermion.
	 */
	const hmc_float mass;

	/**
	 * The system we are operating on.
	 */
	const hardware::System& system;
};


/**
 * A generic fermion matrix (with even-odd preconditioning)
 */
class Fermionmatrix_stagg_eo : public Fermionmatrix_stagg_basic {
public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in) const = 0;

protected:
	Fermionmatrix_stagg_eo(bool herm, hmc_float _mass, const hardware::System& system) : Fermionmatrix_stagg_basic(system, herm, _mass) { };
};


/**
 * Actual fermion matrices (using even-odd)
 */
class MdagM : public Fermionmatrix_stagg_eo {
public:
	MdagM(hmc_float _mass, const hardware::System& system) : Fermionmatrix_stagg_eo(true, _mass, system), tmp(system), tmp2(system) { };
	void operator() (const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in) const override;
	cl_ulong get_flops() const override;
private:
	physics::lattices::Staggeredfield_eo tmp;
	physics::lattices::Staggeredfield_eo tmp2;
};


#if 0
//This is the Wilson code that could be needed in some simulations for staggered.
//In such a case, uncomment and adapt...

 class Fermionmatrix : public Fermionmatrix_basic {
public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const = 0;

protected:
	Fermionmatrix(bool herm, hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_basic(system, herm, _kappa, _mubar) { };

};

/**
 * Actual fermion matrices (no even-odd)
 */
class M : public Fermionmatrix {
public:
	M(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix(false, _kappa, _mubar, system) {  };
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const override;
	cl_ulong get_flops() const override;
};
#endif


}
}
#endif /* _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_STAGG_ */
