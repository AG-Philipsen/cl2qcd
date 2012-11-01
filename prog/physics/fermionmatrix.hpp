#ifndef _FERMIONMATRIXH_
#define _FERMIONMATRIXH_

#include "../opencl_module.h"

#include "../hardware/buffers/plain.hpp"
#include "../hardware/buffers/su3.hpp"
#include "../hardware/buffers/spinor.hpp"
#include "../host_use_timer.h"


/**
 * this is the definition of the class "Fermionmatrix"
 */
namespace physics {

  namespace fermionmatrix{
    /**
     * A generic fermion matrix
     */
    class Fermionmatrix_basic {
    protected:
      Opencl_Module * that;

      Fermionmatrix_basic(Opencl_Module * that) : that(that) { };
      
    private:
      bool is_hermitian;

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
    class Fermionmatrix : public Fermionmatrix_basic{
    protected:
      Fermionmatrix(Opencl_Module * that) : Fermionmatrix_basic(that) { };

    public:
      /**
       * Invoke the matrix function.
       */
      virtual void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;
      
    };
    /**
     * A generic fermion matrix (with even-odd preconditioning)
     */
    class Fermionmatrix_eo : public Fermionmatrix_basic{
    protected:
      Fermionmatrix_eo(Opencl_Module * that) : Fermionmatrix_basic(that) { };
      
    public:
      /**
       * Invoke the matrix function.
       */
      virtual void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;
    };
    
    /**
     * Actual fermion matrices (no even-odd)
     */
    class M : public Fermionmatrix {
    public:
      M(Opencl_Module * that) : Fermionmatrix(that) { };
      void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    class Qplus : public Fermionmatrix {
    public:
      Qplus(Opencl_Module * that) : Fermionmatrix(that) { };
      void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    class Qminus : public Fermionmatrix {
    public:
      Qminus(Opencl_Module * that) : Fermionmatrix(that) { };
      void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    class QplusQminus : public Fermionmatrix {
    public:
      QplusQminus(Opencl_Module * that) : Fermionmatrix(that) { };
      void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    /**
     * Actual fermion matrices (using even-odd)
     */
    class Aee : public Fermionmatrix_eo {
    public:
	Aee(Opencl_Module * that) : Fermionmatrix_eo(that) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
    };
    class Qplus_eo : public Fermionmatrix_eo {
    public:
      Qplus_eo(Opencl_Module * that) : Fermionmatrix_eo(that) { };
      void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    class Qminus_eo : public Fermionmatrix_eo {
    public:
      Qminus_eo(Opencl_Module * that) : Fermionmatrix_eo(that) { };
      void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
    class QplusQminus_eo : public Fermionmatrix_eo {
    public:
      QplusQminus_eo(Opencl_Module * that) : Fermionmatrix_eo(that) { };
      void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
      cl_ulong get_Flops() const override;
      cl_ulong get_Bytes() const override;
    };
  }
}
#endif /* _FERMIONMATRIXH_ */
