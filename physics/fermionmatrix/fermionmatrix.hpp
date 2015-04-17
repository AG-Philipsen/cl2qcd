/** @file
 * Definition of fermionmatrix operations.
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_
#define _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_

#include "../../hardware/code/fermions.hpp"

#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/gaugefield.hpp"

#include "../../hardware/device.hpp"
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

/**
 * A generic fermion matrix
 */
class Fermionmatrix_basic {
public:
	/**
	 * Get if the matrix is hermitian
	 */
	bool is_hermitian() const noexcept;
	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_flops() const = 0;
	/**
	 * Get the net read-write-size used by this function.
	 */
	virtual cl_ulong get_read_write_size() const = 0;

protected:
	Fermionmatrix_basic(const hardware::System& system, bool herm, hmc_float _kappa = ARG_DEF, hmc_float _mubar = ARG_DEF) : _is_hermitian(herm), kappa(_kappa), mubar(_mubar), system(system) { };

	hmc_float get_kappa() const noexcept;
	hmc_float get_mubar() const noexcept;

	/**
	 * Get the system to operate on.
	 */
	const hardware::System& get_system() const noexcept;

private:
	/**
	 * Shows if matrix is hermitian
	 */
	const bool _is_hermitian;

	/*
	 * parameters kappa and mubar
	 */
	const hmc_float kappa;
	const hmc_float mubar;

	/**
	 * The system we are operating on.
	 */
	const hardware::System& system;
};
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
 * A generic fermion matrix (with even-odd preconditioning)
 */
class Fermionmatrix_eo : public Fermionmatrix_basic {
public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const = 0;

protected:
	Fermionmatrix_eo(bool herm, hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_basic(system, herm, _kappa, _mubar) { };
};

/**
 * Actual fermion matrices (no even-odd)
 */
class M : public Fermionmatrix {
public:
	M(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix(false, _kappa, _mubar, system) {  };
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class Qplus : public Fermionmatrix {
public:
	Qplus(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix(false, _kappa, _mubar, system) { };
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class Qminus : public Fermionmatrix {
public:
	Qminus(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix(false, _kappa, _mubar, system) { };
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class QplusQminus : public Fermionmatrix {
public:
	QplusQminus(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix(true, _kappa, _mubar, system), q_plus(_kappa, _mubar, system), q_minus(_kappa, _mubar, system), tmp(system) { };
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Qplus q_plus;
	const Qminus q_minus;
	const physics::lattices::Spinorfield tmp;
};
/**
 * Actual fermion matrices (using even-odd)
 */
class Aee : public Fermionmatrix_eo {
public:
	Aee(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), tmp(system), tmp2(system) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_AND_gamma5_eo : public Fermionmatrix_eo {
public:
	Aee_AND_gamma5_eo(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), tmp(system), tmp2(system) { };
	void operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_minus : public Fermionmatrix_eo {
public:
	Aee_minus(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), tmp(system), tmp2(system) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_minus_AND_gamma5_eo : public Fermionmatrix_eo {
public:
	Aee_minus_AND_gamma5_eo(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), tmp(system), tmp2(system) { };
	void operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Qplus_eo : public Fermionmatrix_eo {
public:
	Qplus_eo(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), aee(_kappa, _mubar, system), aee_AND_gamma5_eo(_kappa, _mubar, system) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Aee aee;
	const Aee_AND_gamma5_eo aee_AND_gamma5_eo;
};
class Qminus_eo : public Fermionmatrix_eo {
public:
	Qminus_eo(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(false, _kappa, _mubar, system), aee_minus(_kappa, _mubar, system), aee_minus_AND_gamma5_eo(_kappa, _mubar, system) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Aee_minus aee_minus;
	const Aee_minus_AND_gamma5_eo aee_minus_AND_gamma5_eo;
};
class QplusQminus_eo : public Fermionmatrix_eo {
public:
	QplusQminus_eo(hmc_float _kappa, hmc_float _mubar, const hardware::System& system) : Fermionmatrix_eo(true, _kappa, _mubar, system), q_plus(_kappa, _mubar, system), q_minus(_kappa, _mubar, system), tmp(system) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Qplus_eo q_plus;
	const Qminus_eo q_minus;
	physics::lattices::Spinorfield_eo tmp;
};
}
}
#endif /* _PHYSICS_FERMIONMATRIX_FERMIONMATRIX_ */
