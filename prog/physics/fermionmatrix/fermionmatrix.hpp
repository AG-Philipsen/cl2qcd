/** @file
 * Definition of fermionmatrix operations.
 */

#ifndef _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_
#define _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_

#include "../hardware/code/fermions.hpp"

#include "../host_use_timer.h"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/gaugefield.hpp"

#include "../hardware/device.hpp"
/**
 * this is the definition of the class "Fermionmatrix"
 */
namespace physics {

namespace fermionmatrix {

/*
 * Explicit Fermion operations
 */
void M_wilson(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF);
void M_tm_plus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void M_tm_minus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void M_tm_inverse_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
void M_tm_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
void M_tm_inverse_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
void M_tm_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar = ARG_DEF);
void dslash(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, int evenodd, hmc_float kappa = ARG_DEF);

/*
 * Compound fermion matrix operations
 */
void M(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void Qplus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void Qminus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void QplusQminus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);

void Qplus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void Qminus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void QplusQminus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void Aee(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);
void Aee_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, hmc_float kappa, hmc_float mubar, const hardware::System& system);

/**
 * A generic fermion matrix
 */
class Fermionmatrix_basic {
protected:
	hardware::code::Fermions * that;

	Fermionmatrix_basic(hardware::code::Fermions * that, bool herm, hmc_float _kappa = ARG_DEF, hmc_float _mubar = ARG_DEF) : that(that), is_hermitian(herm), kappa(_kappa), mubar(_mubar) { };

	/**
	 * Shows if matrix is hermitian
	 */
	const bool is_hermitian;

	/**
	 * parameters kappa and mubar
	 */
	const hmc_float kappa;
	const hmc_float mubar;

public:
	/**
	 * Get if the matrix is hermitian
	 */
	bool get_is_hermitian() const noexcept;
	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_Flops() const = 0;

	/**
	 * Get the net bytes read / written by this function.
	 */
	virtual cl_ulong get_Bytes() const = 0;
};
class Fermionmatrix : public Fermionmatrix_basic {
protected:
	Fermionmatrix(hardware::code::Fermions * that, bool herm, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_basic(that, herm, _kappa, _mubar) { };

public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const = 0;

};
/**
 * A generic fermion matrix (with even-odd preconditioning)
 */
class Fermionmatrix_eo : public Fermionmatrix_basic {
protected:
	Fermionmatrix_eo(hardware::code::Fermions * that, bool herm, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_basic(that, herm, _kappa, _mubar) { };

public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const = 0;
};

/**
 * Actual fermion matrices (no even-odd)
 */
class M : public Fermionmatrix {
public:
	M(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) {  };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus : public Fermionmatrix {
public:
	Qplus(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus : public Fermionmatrix {
public:
	Qminus(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus : public Fermionmatrix {
public:
	QplusQminus(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, true, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
/**
 * Actual fermion matrices (using even-odd)
 */
class Aee : public Fermionmatrix_eo {
public:
	Aee(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus_eo : public Fermionmatrix_eo {
public:
	Qplus_eo(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus_eo : public Fermionmatrix_eo {
public:
	Qminus_eo(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus_eo : public Fermionmatrix_eo {
public:
	QplusQminus_eo(hardware::code::Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, true, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
}
}
#endif /* _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_ */
