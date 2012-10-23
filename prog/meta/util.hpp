/** @file
 * Generic utility functions
 */

#ifndef _META_UTIL_
#define _META_UTIL_

#include "inputparameters.hpp"

#include "../types.h"
#include <iostream>

namespace meta {
	size_t get_volspace(const Inputparameters&);
	size_t get_vol4d(const Inputparameters&);
	bool get_use_rectangles(const Inputparameters& params);
	hmc_float get_mubar(const Inputparameters& params);
	hmc_float get_mubar_mp(const Inputparameters& params);
	/**
	 * @deprecated should be done by the spinorfield class
	 */
	size_t get_spinorfieldsize(const Inputparameters& params);
	/**
	 * @deprecated should be done by the spinorfield class
	 */
	size_t get_eoprec_spinorfieldsize(const Inputparameters& params);
	size_t get_float_size(const Inputparameters& params);
	size_t get_mat_size(const Inputparameters& params);
	size_t get_plaq_norm(const Inputparameters& params);
	size_t get_tplaq_norm(const Inputparameters& params);
	size_t get_splaq_norm(const Inputparameters& params);
	size_t get_rect_norm(const Inputparameters& params);
	size_t get_poly_norm(const Inputparameters& params);
	size_t get_flop_complex_mult() noexcept;
	size_t get_flop_su3_su3() noexcept;
	size_t get_flop_su3trace() noexcept;
	size_t get_flop_su3_su3vec() noexcept;
	size_t get_su3algebrasize() noexcept;
	double get_c0(const Inputparameters& params);
	double get_c1(const Inputparameters& params);
	double get_xi_0(const Inputparameters& params);
	size_t get_flop_spinor_spinor() noexcept;
	size_t get_flop_spinor_sqnorm() noexcept;
	void print_info_gaugeobservables(const char* progname, const Inputparameters& params);
	void print_info_gaugeobservables(const char* progname, std::ostream* os, const Inputparameters& params);
	void print_info_heatbath(const char* progname, const Inputparameters& params);
	void print_info_heatbath(const char* progname, std::ostream* os, const Inputparameters& params);
	void print_info_tkkappa(const char* progname, std::ostream* os, const Inputparameters& params);
	void print_info_tkkappa(const char* progname, const Inputparameters& params);
	void print_info_inverter(const char* progname, const Inputparameters& params);
	void print_info_inverter(const char* progname, std::ostream* os, const Inputparameters& params);
	void print_info_hmc(const char* progname, const Inputparameters& params);
	void print_info_hmc(const char* progname, std::ostream* os, const Inputparameters& params);
  std::string get_ferm_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
  std::string get_gauge_obs_file_name(const Inputparameters& parameters, std::string conf_name) noexcept;
}

#endif /* META_UTIL_ */
