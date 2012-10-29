/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEFERMIONSH_
#define _OPENCLMODULEFERMIONSH_

#include "opencl_module.h"

#include "hardware/buffers/plain.hpp"
#include "hardware/buffers/su3.hpp"
#include "hardware/buffers/spinor.hpp"
#include "host_use_timer.h"

class Opencl_Module_Fermions;

/**
 * this is a workaround to be able to pass a (fermionmatrix-)function, which are methods in this class,
 * to another function inside this class.
 * This type points to a helper-function, which then calls the wanted function.
 */
class Matrix_Function {
protected:
	Opencl_Module_Fermions * that;

	Matrix_Function(Opencl_Module_Fermions * that) : that(that) { };

public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;

	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_Flops() const = 0;

	/**
	 * Get the net bytes read / written by this function.
	 */
	virtual cl_ulong get_Bytes() const = 0;
};

class M : public Matrix_Function {
public:
	M(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus : public Matrix_Function {
public:
	Qplus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus : public Matrix_Function {
public:
	Qminus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus : public Matrix_Function {
public:
	QplusQminus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
/**
 * this is a workaround to be able to pass a (fermionmatrix-)function, which are methods in this class,
 * to another function inside this class.
 * This type points to a helper-function, which then calls the wanted function.
 */
class Matrix_Function_eo {
protected:
	Opencl_Module_Fermions * that;

	Matrix_Function_eo(Opencl_Module_Fermions * that) : that(that) { };

public:
	/**
	 * Invoke the matrix function.
	 */
	virtual void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;

	/**
	 * Get the net flops performed by this function.
	 */
	virtual cl_ulong get_Flops() const = 0;

	/**
	 * Get the net bytes read / written by this function.
	 */
	virtual cl_ulong get_Bytes() const = 0;
};

class Aee : public Matrix_Function_eo {
public:
	Aee(Opencl_Module_Fermions * that) : Matrix_Function_eo(that) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qplus_eo : public Matrix_Function_eo {
public:
	Qplus_eo(Opencl_Module_Fermions * that) : Matrix_Function_eo(that) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class Qminus_eo : public Matrix_Function_eo {
public:
	Qminus_eo(Opencl_Module_Fermions * that) : Matrix_Function_eo(that) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};
class QplusQminus_eo : public Matrix_Function_eo {
public:
	QplusQminus_eo(Opencl_Module_Fermions * that) : Matrix_Function_eo(that) { };
	void operator() (const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const override;
	cl_ulong get_Flops() const override;
	cl_ulong get_Bytes() const override;
};

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Fermions : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Fermions();


	//    fermionmatrix operations
	//    non-eo
	//        compound
	void M(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qplus(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qminus(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void QplusQminus(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//        explicit
	void M_wilson_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF);
	void M_tm_plus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void M_tm_minus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void gamma5_device(const hardware::buffers::Plain<spinor> * inout);
	//    eo
	//        compound
	void Qplus_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qminus_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void QplusQminus_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee_minus(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//        explicit
	void gamma5_eo_device(const hardware::buffers::Spinor * inout);
	void M_tm_inverse_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);
	void M_tm_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);
	void M_tm_inverse_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);
	void M_tm_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);
	void dslash_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF);
	//        merged
	void Aee_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee_minus_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF);
	void dslash_AND_M_tm_inverse_sitediagonal_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
        void M_tm_sitediagonal_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);
        void M_tm_sitediagonal_minus_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF);

	//    solver operations
	//    non-eo
	/// this calls the solver according to parameter settings using the fermionmatrix f
	void solver(const Matrix_Function & f, const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, usetimer * solvertimer);
	/**
	* the solvers return the number of iterations needed if it converged,
	* -1 if it did not converge within cgmax
	* -iter if the algorithm got stuck at some point
	*/
	/// this executes the bicgstab on the device, using the fermionmatrix f
	int bicgstab(const Matrix_Function & f, const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	/// this executes the cg on the device, using the fermionmatrix f
	int cg(const Matrix_Function & f, const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//    eo
	/// this calls the solver according to parameter settings using the fermionmatrix f
	void solver(const Matrix_Function_eo & f, const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, usetimer * solvertimer);
	/// this executes the eo bicgstab on the device, using the fermionmatrix f
	int bicgstab_eo(const Matrix_Function_eo & f, const hardware::buffers::Spinor * inout, const hardware::buffers::Spinor * source, const hardware::buffers::SU3 * gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	int cg_eo(const Matrix_Function_eo & f, const hardware::buffers::Spinor * inout, const hardware::buffers::Spinor * source, const hardware::buffers::SU3 * gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

	/////////////////////////////////////////////////
	//functions to get private variables
	const hardware::buffers::Plain<spinor> * get_inout();
	const hardware::buffers::Plain<spinor> * get_source();
	const hardware::buffers::Plain<spinor> * get_tmp();

	const hardware::buffers::Spinor * get_inout_eo();
	const hardware::buffers::Spinor * get_tmp_eo_1();
	const hardware::buffers::Spinor * get_tmp_eo_2();
	const hardware::buffers::Spinor * get_source_even();
	const hardware::buffers::Spinor * get_source_odd();

	const hardware::buffers::Plain<hmc_complex> * get_clmem_minusone();
	const hardware::buffers::Plain<hmc_complex> * get_clmem_one();

	/**
	 * This is used to print the squarenorm of an inverter-solution in debug-mode
	 * @param eo has to be set to true/false when evenodd is used/not used.
	 * @param msg message to be printed right before before the actual squarenorm is printed
	 */
	hmc_float print_info_inv_field(const hardware::buffers::Plain<spinor> * in, bool eo, std::string msg);

	/**
	 * This is used to print the squarenorm of an inverter-solution in debug-mode
	 * @param eo has to be set to true/false when evenodd is used/not used.
	 * @param msg message to be printed right before before the actual squarenorm is printed
	 */
	hmc_float print_info_inv_field(const hardware::buffers::Spinor * in, bool eo, std::string msg);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;

	ClSourcePackage get_sources() const noexcept;

protected:
	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

private:
	/**
	 * Default constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Fermions(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	void clear_kernels();

	////////////////////////////////////
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M_wilson;
	cl_kernel gamma5;
	cl_kernel M_tm_plus;
	cl_kernel M_tm_minus;
	cl_kernel gamma5_eo;
	cl_kernel M_tm_sitediagonal;
	cl_kernel M_tm_inverse_sitediagonal;
	cl_kernel M_tm_sitediagonal_minus;
	cl_kernel M_tm_inverse_sitediagonal_minus;
	cl_kernel dslash_eo;
	cl_kernel dslash_AND_gamma5_eo;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_eo;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_minus_eo;
	cl_kernel M_tm_sitediagonal_AND_gamma5_eo;
	cl_kernel M_tm_sitediagonal_minus_AND_gamma5_eo;
	//CP: variables for normal solver
	const hardware::buffers::Plain<spinor> clmem_inout;
	const hardware::buffers::Plain<spinor> clmem_source;
	const hardware::buffers::Plain<spinor> clmem_rn;
	const hardware::buffers::Plain<spinor> clmem_rhat;
	const hardware::buffers::Plain<spinor> clmem_v;
	const hardware::buffers::Plain<spinor> clmem_p;
	const hardware::buffers::Plain<spinor> clmem_s;
	const hardware::buffers::Plain<spinor> clmem_t;
	const hardware::buffers::Plain<spinor> clmem_aux;
	//this is needed in QplusQminus as a temporary field
	const hardware::buffers::Plain<spinor> clmem_tmp;

	//CP: variables for eo solver
	const hardware::buffers::Spinor clmem_inout_eo;
	const hardware::buffers::Spinor clmem_source_even;
	const hardware::buffers::Spinor clmem_source_odd;
	const hardware::buffers::Spinor clmem_rn_eo;
	const hardware::buffers::Spinor clmem_rhat_eo;
	const hardware::buffers::Spinor clmem_v_eo;
	const hardware::buffers::Spinor clmem_p_eo;
	const hardware::buffers::Spinor clmem_s_eo;
	const hardware::buffers::Spinor clmem_t_eo;
	const hardware::buffers::Spinor clmem_aux_eo;
	const hardware::buffers::Spinor clmem_tmp_eo_1;
	const hardware::buffers::Spinor clmem_tmp_eo_2;

	const hardware::buffers::Plain<hmc_complex> clmem_rho;
	const hardware::buffers::Plain<hmc_complex> clmem_rho_next;
	const hardware::buffers::Plain<hmc_complex> clmem_alpha;
	const hardware::buffers::Plain<hmc_complex> clmem_omega;
	const hardware::buffers::Plain<hmc_complex> clmem_beta;
	const hardware::buffers::Plain<hmc_complex> clmem_tmp1;
	const hardware::buffers::Plain<hmc_complex> clmem_tmp2;
	const hardware::buffers::Plain<hmc_complex> clmem_one;
	const hardware::buffers::Plain<hmc_complex> clmem_minusone;

	ClSourcePackage sources;
};

#endif //OPENCLMODULEFERMIONSH
