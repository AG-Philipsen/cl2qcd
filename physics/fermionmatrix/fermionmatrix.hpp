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
#include "../lattices/latticesInterfaces.hpp"

/**
 * this is the definition of the class "Fermionmatrix"
 */
namespace physics {

namespace fermionmatrix {

/*
 * Explicit Fermion operations
 */
void M_wilson(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa);
void M_tm_plus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar);
void M_tm_minus(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, hmc_float kappa, hmc_float mubar);
void M_tm_inverse_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar);
void M_tm_sitediagonal(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar);
void M_tm_inverse_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar);
void M_tm_sitediagonal_minus(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Spinorfield_eo& in, hmc_float mubar);
void dslash(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, int evenodd, hmc_float kappa);

/**
 * A generic fermion matrix
 */
class Fermionmatrix_basic {
public:
	bool isHermitian() const noexcept;
	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_flops() const = 0;
	/**
	 * Get the net read-write-size used by this function.
	 */
	virtual cl_ulong get_read_write_size() const = 0;
	virtual ~Fermionmatrix_basic() {};

protected:
	Fermionmatrix_basic(const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface,
	                    bool herm)
        : fermionmatrixParametersInterface(fermionmatrixParametersInterface), _is_hermitian(herm), system(system) {};

	/**
	 * Get the system to operate on.
	 */
	const hardware::System& get_system() const noexcept;
	const FermionmatrixParametersInterface& fermionmatrixParametersInterface;

private:
	/**
	 * Shows if matrix is hermitian
	 */
	const bool _is_hermitian;

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
	virtual void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& in, const physics::AdditionalParameters& additionalParameters) const = 0;
	virtual ~Fermionmatrix() {};

protected:
	Fermionmatrix(bool herm, const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface)
        : Fermionmatrix_basic(system, fermionmatrixParametersInterface, herm) {};

};
/**
 * A generic fermion matrix (with even-odd preconditioning)
 */
class Fermionmatrix_eo : public Fermionmatrix_basic {
public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const = 0;
	virtual ~Fermionmatrix_eo() {};

protected:
	Fermionmatrix_eo(bool herm, const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface)
        : Fermionmatrix_basic(system, fermionmatrixParametersInterface, herm) {};
};

/**
 * Actual fermion matrices (no even-odd)
 */
class M final : public Fermionmatrix {
public:
	M(const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface)
        : Fermionmatrix(false, system, fermionmatrixParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class Qplus final : public Fermionmatrix {
public:
	Qplus(const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface)
        : Fermionmatrix(false, system, fermionmatrixParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class Qminus final : public Fermionmatrix {
public:
	Qminus(const hardware::System& system, const FermionmatrixParametersInterface& fermionmatrixParametersInterface)
        : Fermionmatrix(false, system, fermionmatrixParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
};
class QplusQminus final : public Fermionmatrix {
public:
	QplusQminus(const hardware::System& system, const FermionParametersInterface& fermionParametersInterface)
        : Fermionmatrix(true,system, fermionParametersInterface), q_plus(system, fermionParametersInterface),
          q_minus(system, fermionParametersInterface), tmp(system, fermionParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield& in, const physics::AdditionalParameters& additionalParameters) const override;
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
class Aee final : public Fermionmatrix_eo {
public:
	Aee(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          tmp(system, fermionEoParametersInterface), tmp2(system, fermionEoParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_AND_gamma5_eo final : public Fermionmatrix_eo {
public:
	Aee_AND_gamma5_eo(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          tmp(system, fermionEoParametersInterface), tmp2(system, fermionEoParametersInterface) {};
	void operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_minus final : public Fermionmatrix_eo {
public:
	Aee_minus(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          tmp(system, fermionEoParametersInterface), tmp2(system, fermionEoParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Aee_minus_AND_gamma5_eo final : public Fermionmatrix_eo {
public:
	Aee_minus_AND_gamma5_eo(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          tmp(system, fermionEoParametersInterface), tmp2(system, fermionEoParametersInterface) {};
	void operator()(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	physics::lattices::Spinorfield_eo tmp;
	physics::lattices::Spinorfield_eo tmp2;
};
class Qplus_eo final : public Fermionmatrix_eo {
public:
	Qplus_eo(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          aee(system, fermionEoParametersInterface),
          aee_AND_gamma5_eo(system, fermionEoParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Aee aee;
	const Aee_AND_gamma5_eo aee_AND_gamma5_eo;
};
class Qminus_eo : public Fermionmatrix_eo {
public:
	Qminus_eo(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(false, system, fermionEoParametersInterface),
          aee_minus(system, fermionEoParametersInterface),
          aee_minus_AND_gamma5_eo(system, fermionEoParametersInterface) {};
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
	cl_ulong get_flops() const override;
	cl_ulong get_read_write_size() const override;
private:
	const Aee_minus aee_minus;
	const Aee_minus_AND_gamma5_eo aee_minus_AND_gamma5_eo;
};
class QplusQminus_eo final : public Fermionmatrix_eo {
public:
	QplusQminus_eo(const hardware::System& system, const FermionEoParametersInterface& fermionEoParametersInterface)
        : Fermionmatrix_eo(true, system, fermionEoParametersInterface),
          q_plus(system, fermionEoParametersInterface),
          q_minus(system, fermionEoParametersInterface),
          tmp(system, fermionEoParametersInterface) { };
	void operator() (const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
	                 const physics::lattices::Spinorfield_eo& in, const physics::AdditionalParameters& additionalParameters) const override;
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
