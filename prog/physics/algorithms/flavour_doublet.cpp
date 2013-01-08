/** @file
 * Implementation of the flavour doublets algorithms
 */

#include "flavour_doublet.hpp"

#include <stdexcept>
#include <cassert>
#include <fstream>
#include "../meta/util.hpp"
#include "../lattices/util.hpp"

static hardware::buffers::Plain<spinor> * merge_spinorfields(const std::vector<const physics::lattices::Spinorfield*>& fields, const size_t device_idx, hardware::Device * device);

static hmc_complex flavour_doublet_chiral_condensate_std(const std::vector<const physics::lattices::Spinorfield*>& solved_fields, const std::vector<const physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system);
static hmc_complex flavour_doublet_chiral_condensate_tm(const std::vector<const physics::lattices::Spinorfield*>& solved_fields, std::string pbp_fn, int number, const hardware::System& system);
static size_t get_num_corr_entries(const meta::Inputparameters& params);

void physics::algorithms::flavour_doublet_correlators(const std::vector<const physics::lattices::Spinorfield*>& result, const std::vector<const physics::lattices::Spinorfield*>& sources, std::ostream& of, const meta::Inputparameters& parameters)
{
	using namespace std;

	auto result_ps = calculate_correlator("ps", result, sources, parameters);
	auto result_sc = calculate_correlator("sc", result, sources, parameters);
	auto result_vx = calculate_correlator("vx", result, sources, parameters);
	auto result_vy = calculate_correlator("vy", result, sources, parameters);
	auto result_vz = calculate_correlator("vz", result, sources, parameters);
	auto result_ax = calculate_correlator("ax", result, sources, parameters);
	auto result_ay = calculate_correlator("ay", result, sources, parameters);
	auto result_az = calculate_correlator("az", result, sources, parameters);

	if(parameters.get_print_to_screen() )
		meta::print_info_flavour_doublet_correlators(parameters);

	meta::print_info_flavour_doublet_correlators(&of, parameters);

	// @todo One could also implement to write all results on screen if wanted
	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	for(int j = 0; j < result_ps.size(); j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << result_ps[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << result_ps[j] << endl;
	}

	//the scalar (J=0, P=0)
	for(int j = 0; j < result_sc.size(); j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << result_sc[j] << endl;
	}

	//the vector (J=1, P=1)
	if(result_vx.size() != result_vy.size() || result_vx.size() != result_vz.size()) {
		throw Print_Error_Message("Internal error: Vector correlators are not of equal length");
	}
	for(int j = 0; j < result_vx.size(); j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (result_vx[j] + result_vy[j] + result_vz[j]) / 3. << "\t" << result_vx[j] << "\t" << result_vy[j] << "\t" << result_vz[j] << endl;
	}

	//the axial vector (J=1, P=0)
	if(result_ax.size() != result_ay.size() || result_ax.size() != result_az.size()) {
		throw Print_Error_Message("Internal error: Vector correlators are not of equal length");
	}
	for(int j = 0; j < result_ax.size(); j++) {
		of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (result_ax[j] + result_ay[j] + result_az[j]) / 3. << "\t" << result_ax[j] << "\t" << result_ay[j] << "\t" << result_az[j] << endl;
	}

	of << endl;
}

void physics::algorithms::flavour_doublet_chiral_condensate(const std::vector<const physics::lattices::Spinorfield*>& inverted, const std::vector<const physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system)
{
	using namespace std;

	auto params = system.get_inputparameters();

	hmc_complex result;
	if(params.get_pbp_version() == meta::Inputparameters::std) {
		result = flavour_doublet_chiral_condensate_std(inverted, sources, pbp_fn, number, system);
	} else if(params.get_fermact() == meta::Inputparameters::twistedmass && params.get_pbp_version() == meta::Inputparameters::tm_one_end_trick ) {
		result = flavour_doublet_chiral_condensate_tm(inverted, pbp_fn, number, system);
	} else {
		throw std::invalid_argument("No valid chiral condensate version has ben selected.");
	}

	// Output
	ofstream of;
	of.open(pbp_fn.c_str(), ios_base::app);
	if(!of.is_open()) {
		throw File_Exception(pbp_fn);
	}
	logger.info() << "chiral condensate:" ;
	logger.info() << number << "\t" << scientific << setprecision(14) << result.re << "\t" << result.im;
	of << number << "\t" << scientific << setprecision(14) << result.re << "\t" << result.im << endl;
}

static hardware::buffers::Plain<spinor> * merge_spinorfields(const std::vector<const physics::lattices::Spinorfield*>& fields, const size_t device_idx, hardware::Device * device)
{
	size_t total_elems = 0;
for(auto field: fields) {
		total_elems += field->get_buffers().at(device_idx)->get_elements();
	}
	hardware::buffers::Plain<spinor> * result = new hardware::buffers::Plain<spinor>(total_elems, device);

	size_t offset = 0;
for(auto field: fields) {
		auto buffer = field->get_buffers().at(device_idx);
		size_t elems = buffer->get_elements();
		result->copyDataBlock(buffer, offset);
		offset += elems;
	}
	return result;
}

static hmc_complex flavour_doublet_chiral_condensate_std(const std::vector<const physics::lattices::Spinorfield*>& solved_fields, const std::vector<const physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system)
{
	using namespace physics::lattices;

	auto params = system.get_inputparameters();

	hmc_complex result = {0., 0.};
	/**
	 * In the pure Wilson case one can evaluate <pbp> with stochastic estimators according to:
	 * <pbp> = <ubu + dbd> = 2<ubu> = 2 Tr_(space, colour, dirac) ( D^-1 )
	 auto spinor_code = solver->get_device()->get_spinor_code(); *       = lim_r->inf 2/r (Xi_r, Phi_r)
	 * where the estimators satisfy
	 * D^-1(x,y)_(a,b, A,B) = lim_r->inf Phi_r(x)_a,A (Xi_r(y)_b,B)^dagger
	 * and Phi fulfills
	 * D Phi = Xi
	 * (X,Y) denotes the normal scalar product
	 * In the twisted-mass case one can evaluate <pbp> with stochastic estimators similar to the pure wilson case.
	 * However, one first has to switch to the twisted basis:
	 * <pbp> -> <chibar i gamma_5 tau_3 chi>
	 *       = <ub i gamma_5 u> - <db i gamma_5 d>
	 *       = Tr( i gamma_5 (D^-1_u - D^-1_d ) )
	 *       = Tr( i gamma_5 (D^-1_u - gamma_5 D^-1_u^dagger gamma_5) )
	 *       = Tr( i gamma_5 (D^-1_u -  D^-1_u^dagger ) )
	 *       = 2 Im Tr ( gamma_5 D^-1_u)
	 *       = lim_r->inf 2/r  (gamma_5 Xi_r, Phi_r)
	 * NOTE: The basic difference compared to the pure Wilson case is only the gamma_5 and that one takes the imaginary part!
	 */
	// Need 2 spinors at once..
	logger.debug() << "init buffers for chiral condensate calculation...";
	Spinorfield phi(system);
	Spinorfield xi(system);
	assert(solved_fields.size() == sources.size());
	for(size_t i = 0; i < solved_fields.size(); ++i) {
		copyData(&phi, solved_fields[i]);
		copyData(&xi, sources[i]);

		if(params.get_fermact() == meta::Inputparameters::twistedmass) {
			xi.gamma5();
		}

		hmc_complex tmp = scalar_product(xi, phi);
		switch(params.get_fermact()) {
			case  meta::Inputparameters::wilson:
				result.re += tmp.re;
				result.im += tmp.im;
				break;
			case meta::Inputparameters::twistedmass:
				//NOTE: Here I just swap the assignments
				result.re += tmp.im;
				result.im += tmp.re;
				break;
			default:
				throw std::invalid_argument("chiral condensate not implemented for given fermion action");
		}
	}
	//Normalization: The 2/r from above plus 1/VOL4D
	hmc_float norm = 2. / params.get_num_sources() / meta::get_vol4d(params);
	result.re *= norm;
	result.im *= norm;

	return result;
}

static hmc_complex flavour_doublet_chiral_condensate_tm(const std::vector<const physics::lattices::Spinorfield*>& solved_fields, std::string pbp_fn, int number, const hardware::System& system)
{
	hmc_complex result = {0., 0.};
	/**
	 * For twisted-mass fermions one can also employ the one-end trick, which origins from
	 * D_d - D_u = - 4 i kappa amu gamma_5 <-> D^-1_u - D^-1_d = - 4 i kappa amu gamma_5 (D^-1_u)^dagger D^-1_u
	 * With this, the chiral condensate is:
	 * <pbp> = ... = Tr( i gamma_5 (D^-1_u - D^-1_d ) )
	 *       = - 4 kappa amu lim_r->inf 1/R (Phi_r, Phi_r)
	 * NOTE: Here one only needs Phi...
	 */
	logger.debug() << "init buffers for chiral condensate calculation...";
for(auto phi: solved_fields) {
		hmc_float tmp = squarenorm(*phi);
		result.re += tmp;
		result.im += 0;
	}
	//Normalization: The 1/r from above plus 1/VOL4D and the additional factor - 4 kappa amu ( = 2 mubar )
	auto params = system.get_inputparameters();
	hmc_float norm = 1. / params.get_num_sources() / meta::get_vol4d(params) * (-2.) * meta::get_mubar(params);
	result.re *= norm;
	result.im *= norm;
	return result;
}

std::vector<hmc_float> physics::algorithms::calculate_correlator(std::string type, const std::vector<const physics::lattices::Spinorfield*>& corr, const std::vector<const physics::lattices::Spinorfield*>& sources, const meta::Inputparameters& params)
{
	// assert single device
	auto first_field_buffers = corr.at(0)->get_buffers();
	// require single device
	assert(first_field_buffers.size() == 1);
	hardware::Device * device = first_field_buffers.at(0)->get_device();
	auto code = device->get_correlator_code();


	const size_t num_corr_entries = get_num_corr_entries(params);
	const hardware::buffers::Plain<hmc_float> result(num_corr_entries, device);
	result.clear();

	// for each source
	if(corr.size() != sources.size()) {
		throw std::invalid_argument("Correlated and source fields need to be of the same size.");
	}
	// TODO adjust correlator kernels!
	if(params.get_sourcetype() == meta::Inputparameters::point) {
		auto merged_corrs = merge_spinorfields(corr, 0, device);
		code->correlator(code->get_correlator_kernel(type), &result, merged_corrs);
		delete merged_corrs;
	} else {
		auto merged_corrs = merge_spinorfields(corr, 0, device);
		auto merged_sources = merge_spinorfields(sources, 0, device);
		code->correlator(code->get_correlator_kernel(type), &result, merged_corrs, merged_sources);
		delete merged_sources;
		delete merged_corrs;
	}

	std::vector<hmc_float> out(num_corr_entries);
	result.dump(out.data());
	return out;
}

static size_t get_num_corr_entries(const meta::Inputparameters& parameters)
{
	switch (parameters.get_corr_dir()) {
		case 0 :
			return parameters.get_ntime();
		case 3 :
			return parameters.get_nspace();
		default :
			std::stringstream errmsg;
			errmsg << "Correlator direction " << parameters.get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}
}
