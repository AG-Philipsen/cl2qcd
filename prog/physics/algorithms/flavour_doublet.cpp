/** @file
 * Implementation of the flavour doublets algorithms
 */

#include "flavour_doublet.hpp"

#include <stdexcept>
#include <cassert>
#include <fstream>
#include <../meta/util.hpp>

static hardware::buffers::Plain<spinor> * merge_spinorfields(const std::vector<const physics::lattices::Spinorfield*>& fields, const size_t device_idx, hardware::Device * device);

void calculate_correlator(hmc_float* out, size_t num_corr_entries, std::string type, const std::vector<const physics::lattices::Spinorfield*>& corr, const std::vector<const physics::lattices::Spinorfield*>& sources, const meta::Inputparameters& params)
{
	// assert single device
	auto first_field_buffers = corr.at(0)->get_buffers();
	// require single device
	assert(first_field_buffers.size() == 1);
	hardware::Device * device = first_field_buffers.at(0)->get_device();
	auto code = device->get_correlator_code();


	const hardware::buffers::Plain<hmc_float> result(num_corr_entries, device);
	result.clear();

	// for each source
	if(corr.size() != sources.size()) {
		throw std::invalid_argument("Correlated and source fields need to be of the same size.");
	}
	// TODO adjust correlator kernels!
	if(params.get_sourcetype() == meta::Inputparameters::point) {
		auto merged_corrs = merge_spinorfields(corr, 0, device);
		code->correlator_device(code->get_correlator_kernel(type), merged_corrs, &result);
		delete merged_corrs;
	} else {
		auto merged_corrs = merge_spinorfields(corr, 0, device);
		auto merged_sources = merge_spinorfields(sources, 0, device);
		code->correlator_device(code->get_correlator_kernel(type), merged_corrs, merged_sources, &result);
		delete merged_sources;
		delete merged_corrs;
	}

	result.dump(out);
}

void physics::algorithms::flavour_doublet_correlators(const std::vector<const physics::lattices::Spinorfield*>& result, const std::vector<const physics::lattices::Spinorfield*>& sources, std::string corr_fn, const meta::Inputparameters& parameters)
{
	using namespace std;

	int num_corr_entries =  0;
	switch (parameters.get_corr_dir()) {
		case 0 :
			num_corr_entries = parameters.get_ntime();
			break;
		case 3 :
			num_corr_entries = parameters.get_nspace();
			break;
		default :
			stringstream errmsg;
			errmsg << "Correlator direction " << parameters.get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}

	hmc_float* host_result_ps = new hmc_float [num_corr_entries];
	hmc_float* host_result_sc = new hmc_float [num_corr_entries];
	hmc_float* host_result_vx = new hmc_float [num_corr_entries];
	hmc_float* host_result_vy = new hmc_float [num_corr_entries];
	hmc_float* host_result_vz = new hmc_float [num_corr_entries];
	hmc_float* host_result_ax = new hmc_float [num_corr_entries];
	hmc_float* host_result_ay = new hmc_float [num_corr_entries];
	hmc_float* host_result_az = new hmc_float [num_corr_entries];

	calculate_correlator(host_result_ps, num_corr_entries, "ps", result, sources, parameters);
	calculate_correlator(host_result_sc, num_corr_entries, "sc", result, sources, parameters);
	calculate_correlator(host_result_vx, num_corr_entries, "vx", result, sources, parameters);
	calculate_correlator(host_result_vy, num_corr_entries, "vy", result, sources, parameters);
	calculate_correlator(host_result_vz, num_corr_entries, "vz", result, sources, parameters);
	calculate_correlator(host_result_ax, num_corr_entries, "ax", result, sources, parameters);
	calculate_correlator(host_result_ay, num_corr_entries, "ay", result, sources, parameters);
	calculate_correlator(host_result_az, num_corr_entries, "az", result, sources, parameters);

	if(parameters.get_print_to_screen() )
		meta::print_info_flavour_doublet_correlators(parameters);

	ofstream of(corr_fn.c_str(), ios_base::app);
	if(of.is_open()) {
		meta::print_info_flavour_doublet_correlators(&of, parameters);
	} else {
		throw File_Exception(corr_fn);
	}

	// @todo One could also implement to write all results on screen if wanted
	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	for(int j = 0; j < num_corr_entries; j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << host_result_ps[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << host_result_ps[j] << endl;
	}

	//the scalar (J=0, P=0)
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << host_result_sc[j] << endl;
	}

	//the vector (J=1, P=1)
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (host_result_vx[j] + host_result_vy[j] + host_result_vz[j]) / 3. << "\t" << host_result_vx[j] << "\t" << host_result_vy[j] << "\t" << host_result_vz[j] << endl;
	}

	//the axial vector (J=1, P=0)
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (host_result_ax[j] + host_result_ay[j] + host_result_az[j]) / 3. << "\t" << host_result_ax[j] << "\t" << host_result_ay[j] << "\t" << host_result_az[j] << endl;
	}

	of << endl;
	of.close();
	delete [] host_result_ps;
	delete [] host_result_sc;
	delete [] host_result_vx;
	delete [] host_result_vy;
	delete [] host_result_vz;
	delete [] host_result_ax;
	delete [] host_result_ay;
	delete [] host_result_az;
}

//void hardware::algorithms::flavour_doublet_chiral_condensate(const physics::lattices::Gaugefield& gaugefield, const std::vector<const physics::lattices::Spinorfield*>& result, std::string pbp_fn, int number);

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
		offset += elems;
		result->copyDataBlock(buffer, offset);
	}
	return result;
}
