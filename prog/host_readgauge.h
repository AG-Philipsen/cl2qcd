/** @file
 * Reading of gauge field from files.
 */

#ifndef _READGAUGEH_
#define _READGAUGEH_

#include "types.h"

/**
 * Parser class for a stored gaugefield.
 *
 * Contains metadata of the parsed gaugefield as members.
 */
class sourcefileparameters {
public:
	sourcefileparameters() {
		set_defaults();
	};
	/**
	 * Read gauge configuration from the given file into the given array.
	 *
	 * @param[in] file      The file to read the gauge configuration from
	 * @param[in] precision The precision expected for the gaugefield.
	 * @param[out] array    The loaded gaugefield
	 */
	void readsourcefile(const char * file, int precision, char ** data);
	void set_defaults();
	void val_assign_source(int * out, int in);
	void val_assign_source(hmc_float * out, hmc_float in);
	int lx_source, ly_source, lz_source, lt_source, prec_source, num_entries_source, flavours_source,
	    trajectorynr_source, time_source, time_solver_source, noiter_source;
	hmc_float plaquettevalue_source, beta_source, kappa_source, mu_source, c2_rec_source, mubar_source, epsilonbar_source,
	          epssq_source, kappa_solver_source, mu_solver_source;
	char field_source[100];
	char hmcversion_source[50];
	char date_source[50];
	char solvertype_source[50];
	char hmcversion_solver_source[50];
	char date_solver_source[50];
};

#endif /* _READGAUGEH_ */
