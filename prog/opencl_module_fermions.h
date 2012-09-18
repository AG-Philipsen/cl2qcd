/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEFERMIONSH_
#define _OPENCLMODULEFERMIONSH_

#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "types_fermions.h"
#include "host_use_timer.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_spinors.h"

#include "exceptions.h"

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
	virtual void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const = 0;

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
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class Qplus : public Matrix_Function {
public:
	Qplus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class Qminus : public Matrix_Function {
public:
	Qminus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class QplusQminus : public Matrix_Function {
public:
	QplusQminus(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class Aee : public Matrix_Function {
public:
	Aee(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class Qplus_eo : public Matrix_Function {
public:
	Qplus_eo(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class Qminus_eo : public Matrix_Function {
public:
	Qminus_eo(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};
class QplusQminus_eo : public Matrix_Function {
public:
	QplusQminus_eo(Opencl_Module_Fermions * that) : Matrix_Function(that) { };
	void operator() (cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	cl_ulong get_Flops() const;
	cl_ulong get_Bytes() const;
};

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Fermions : public Opencl_Module_Spinors {
public:
	/**
	 * Empty constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Fermions(const meta::Inputparameters& params)
		: Opencl_Module_Spinors(params) { }


	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(std::stringstream* collect_options);
	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers();
	/**
	 * Collect the buffers related to the solver.
	 */
	void fill_solver_buffers();
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();
	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels();
	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers();
	/**
	 * Clear out the buffers related to the solver
	 */
	virtual void clear_solver_buffers();

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);


	//    fermionmatrix operations
	//    non-eo
	//        compound
	void M(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qplus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qminus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void QplusQminus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//        explicit
	void M_wilson_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF);
	void M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void gamma5_device(cl_mem inout);
	//    eo
	//        compound
	void Qplus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Qminus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void QplusQminus_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee_minus(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//        explicit
	void gamma5_eo_device(cl_mem inout);
	void M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);
	void M_tm_sitediagonal_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);
	void M_tm_inverse_sitediagonal_minus_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);
	void M_tm_sitediagonal_minus_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);
	void dslash_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa = ARG_DEF);
	//        merged
	void Aee_AND_gamma5_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void Aee_minus_AND_gamma5_eo(cl_mem in, cl_mem out, cl_mem gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_gamma5_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa = ARG_DEF);
	void dslash_AND_M_tm_inverse_sitediagonal_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
        void M_tm_sitediagonal_AND_gamma5_eo_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);
        void M_tm_sitediagonal_minus_AND_gamma5_eo_device(cl_mem in, cl_mem out, hmc_float mubar = ARG_DEF);

	//    solver operations
	//    non-eo
	/// this calls the solver according to parameter settings using the fermionmatrix f
	void solver(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, usetimer * solvertimer);
	/**
	* the solvers return the number of iterations needed if it converged,
	* -1 if it did not converge within cgmax
	* -iter if the algorithm got stuck at some point
	*/
	/// this executes the bicgstab on the device, using the fermionmatrix f
	int bicgstab(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	/// this executes the cg on the device, using the fermionmatrix f
	int cg(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	//    eo
	/// this executes the eo bicgstab on the device, using the fermionmatrix f
	int bicgstab_eo(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	int cg_eo(const Matrix_Function & f, cl_mem inout, cl_mem source, cl_mem gf, hmc_float prec, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_inout();
	cl_mem get_clmem_source();
	cl_mem get_clmem_tmp();

	cl_mem get_clmem_inout_eo();
	cl_mem get_clmem_tmp_eo_1();
	cl_mem get_clmem_tmp_eo_2();
	cl_mem get_clmem_source_even();
	cl_mem get_clmem_source_odd();

	cl_mem get_clmem_minusone();
	cl_mem get_clmem_one();

	/**
	 * This is used to print the squarenorm of an inverter-solution in debug-mode
	 * @param eo has to be set to true/false when evenodd is used/not used.
	 * @param msg message to be printed right before before the actual squarenorm is printed
	 */
	hmc_float print_info_inv_field(cl_mem in, bool eo, std::string msg);

	//protected:

#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	//fermionmatrix
	usetimer timer_M_wilson;
	usetimer timer_gamma5;
	usetimer timer_M_tm_plus;
	usetimer timer_M_tm_minus;
	usetimer timer_gamma5_eo;
	usetimer timer_M_tm_sitediagonal;
	usetimer timer_M_tm_inverse_sitediagonal;
	usetimer timer_dslash_eo;
	usetimer timer_M_tm_sitediagonal_minus;
	usetimer timer_M_tm_inverse_sitediagonal_minus;
	usetimer timer_dslash_AND_gamma5_eo;
	usetimer timer_dslash_AND_M_tm_inverse_sitediagonal_eo;
	usetimer timer_dslash_AND_M_tm_inverse_sitediagonal_minus_eo;
	usetimer timer_M_tm_sitediagonal_AND_gamma5_eo;
	usetimer timer_M_tm_sitediagonal_minus_AND_gamma5_eo;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(const char * in);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(std::string filename, int number);

#endif

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const char * in);

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const char * in);

private:
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
	cl_mem clmem_inout;
	cl_mem clmem_rn;
	cl_mem clmem_rhat;
	cl_mem clmem_v;
	cl_mem clmem_p;
	cl_mem clmem_s;
	cl_mem clmem_t;
	cl_mem clmem_aux;
	//this is needed in QplusQminus as a temporary field
	cl_mem clmem_tmp;

	//CP: variables for eo solver
	cl_mem clmem_inout_eo;
	cl_mem clmem_source;
	cl_mem clmem_source_even;
	cl_mem clmem_source_odd;
	cl_mem clmem_rn_eo;
	cl_mem clmem_rhat_eo;
	cl_mem clmem_v_eo;
	cl_mem clmem_p_eo;
	cl_mem clmem_s_eo;
	cl_mem clmem_t_eo;
	cl_mem clmem_aux_eo;
	cl_mem clmem_tmp_eo_1;
	cl_mem clmem_tmp_eo_2;

	cl_mem clmem_rho;
	cl_mem clmem_rho_next;
	cl_mem clmem_alpha;
	cl_mem clmem_omega;
	cl_mem clmem_beta;
	cl_mem clmem_tmp1;
	cl_mem clmem_tmp2;
	cl_mem clmem_one;
	cl_mem clmem_minusone;

	cl_mem clmem_resid;
	cl_mem clmem_trueresid;
};

#endif //OPENCLMODULEFERMIONSH
