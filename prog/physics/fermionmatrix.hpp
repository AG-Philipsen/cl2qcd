#ifndef _FERMIONMATRIXH_
#define _FERMIONMATRIXH_

#include "../hardware/code/fermions.hpp"

#include "../hardware/buffers/plain.hpp"
#include "../hardware/buffers/su3.hpp"
#include "../hardware/buffers/spinor.hpp"
#include "../host_use_timer.h"

#include "../hardware/device.hpp"
/**
 * this is the definition of the class "Fermionmatrix"
 */
namespace physics {

namespace fermionmatrix {
/**
 * A generic fermion matrix
 */
class Fermionmatrix_basic {
protected:
	Opencl_Module_Fermions * that;

	Fermionmatrix_basic(Opencl_Module_Fermions * that, bool herm, hmc_float _kappa = ARG_DEF, hmc_float _mubar = ARG_DEF) : that(that), is_hermitian(herm), kappa(_kappa), mubar(_mubar) { };

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
	Fermionmatrix(Opencl_Module_Fermions * that, bool herm, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_basic(that, herm, _kappa, _mubar) { };

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
	Fermionmatrix_eo(Opencl_Module_Fermions * that, bool herm, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_basic(that, herm, _kappa, _mubar) { };

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
	M(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) {  };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus : public Fermionmatrix {
public:
	Qplus(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus : public Fermionmatrix {
public:
	Qminus(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus : public Fermionmatrix {
public:
	QplusQminus(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix(that, true, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
/**
 * Actual fermion matrices (using even-odd)
 */
class Aee : public Fermionmatrix_eo {
public:
	Aee(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus_eo : public Fermionmatrix_eo {
public:
	Qplus_eo(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus_eo : public Fermionmatrix_eo {
public:
	Qminus_eo(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, false, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus_eo : public Fermionmatrix_eo {
public:
	QplusQminus_eo(Opencl_Module_Fermions * that, hmc_float _kappa, hmc_float _mubar) : Fermionmatrix_eo(that, true, _kappa, _mubar) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
}
}
#endif /* _FERMIONMATRIXH_ */
