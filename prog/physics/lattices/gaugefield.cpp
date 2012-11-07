/** @file
 * Implementation of the physics::lattices::Gaugefield class
 */

#include "gaugefield.hpp"
#include "../meta/version.hpp"
#include "../../meta/util.hpp"
#include "../../logger.hpp"
#include "../../host_operations_gaugefield.h"
#include "../../host_writegaugefield.h"
#include "../../host_readgauge.h"
#include <cassert>

/**
 * Version number.
 *
 * @deprecated move this into some specific header or so
 */
extern std::string const version;

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system);

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, physics::PRNG& prng);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);
static void set_cold(Matrixsu3 * field, size_t elems);
static void set_hot(Matrixsu3 * field, physics::PRNG& prng, size_t elems);
static void copy_gaugefield_to_ildg_format(hmc_float * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters);
static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, hmc_float * gaugefield_tmp, int check, const meta::Inputparameters& parameters);
static void check_sourcefileparameters(const meta::Inputparameters& parameters, const hmc_float, sourcefileparameters& parameters_source);

physics::lattices::Gaugefield::Gaugefield(hardware::System& system, physics::PRNG& prng)
	: system(system), prng(prng), buffers(allocate_buffers(system))
{
	auto parameters = system.get_inputparameters();
	switch(parameters.get_startcondition()) {
		case meta::Inputparameters::start_from_source:
			fill_from_ildg(parameters.get_sourcefile());
			break;
		case meta::Inputparameters::cold_start:
			set_cold(buffers);
			break;
		case meta::Inputparameters::hot_start:
			set_hot(buffers, prng);
			break;
	}
}

physics::lattices::Gaugefield::Gaugefield(hardware::System& system, physics::PRNG& prng, bool hot)
	: system(system), prng(prng), buffers(allocate_buffers(system))
{
	if(hot) {
		set_hot(buffers, prng);
	} else {
		set_cold(buffers);
	}
}

physics::lattices::Gaugefield::Gaugefield(hardware::System& system, physics::PRNG& prng, std::string ildgfile)
	: system(system), prng(prng), buffers(allocate_buffers(system))
{
	fill_from_ildg(ildgfile);
}

void physics::lattices::Gaugefield::fill_from_ildg(std::string ildgfile)
{
	assert(buffers.size() == 1);

	auto parameters = system.get_inputparameters();
	hmc_float * gf_ildg = new hmc_float[NDIM * NC * NC * parameters.get_ntime() * meta::get_volspace(parameters)];
	Matrixsu3 * gf_host = new Matrixsu3[buffers[0]->get_elements()];

	sourcefileparameters parameters_source;
	parameters_source.readsourcefile(ildgfile.c_str(), parameters.get_precision(), &gf_ildg);
	copy_gaugefield_from_ildg_format(gf_host, gf_ildg, parameters_source.num_entries_source, parameters);

	auto device = buffers[0]->get_device();
	device->get_gaugefield_code()->importGaugefield(buffers[0], gf_host);
	device->synchronize();

	delete[] gf_ildg;
	delete[] gf_host;

	hmc_float plaq = plaquette();
	check_sourcefileparameters(parameters, plaq, parameters_source);
}

static  std::vector<const hardware::buffers::SU3 *> allocate_buffers(hardware::System& system)
{
	using hardware::buffers::SU3;

	// only use device 0 for now
	hardware::Device * device = system.get_devices().at(0);
	std::vector<const SU3 *> buffers;
	buffers.push_back(new SU3(meta::get_vol4d(system.get_inputparameters()) * 4, device));
	return buffers;
}

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, physics::PRNG& prng)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_hot(tmp, prng, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_cold(tmp, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

void set_cold(Matrixsu3 * field, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void set_hot(Matrixsu3 * field, physics::PRNG& prng, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = random_matrixsu3();
	}
}

void physics::lattices::Gaugefield::save(int number)
{
	std::string outputfile = meta::create_configuration_name(system.get_inputparameters(), number);
	save(outputfile);
}

void physics::lattices::Gaugefield::save(std::string outputfile)
{
	assert(buffers.size() == 1);

	auto parameters = system.get_inputparameters();
	const size_t NTIME = parameters.get_ntime();
	const size_t gaugefield_buf_size = 2 * NC * NC * NDIM * meta::get_volspace(parameters) * NTIME;
	hmc_float * gaugefield_buf = new hmc_float[gaugefield_buf_size];

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;

	{
		auto dev_buf = buffers[0];
		Matrixsu3 * host_buf = new Matrixsu3[dev_buf->get_elements()];
		auto device = dev_buf->get_device();
		device->get_gaugefield_code()->exportGaugefield(host_buf, dev_buf);
		device->synchronize();

		copy_gaugefield_to_ildg_format(gaugefield_buf, host_buf, parameters);

		delete host_buf;
	}

	hmc_float plaq = plaquette();

	int number = 0;

	const size_t NSPACE = parameters.get_nspace();

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NTIME, parameters.get_precision(), number, plaq, parameters.get_beta(), parameters.get_kappa(), parameters.get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	delete[] gaugefield_buf;
}

static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, hmc_float * gaugefield_tmp, int check, const meta::Inputparameters& parameters)
{
	//little check if arrays are big enough
	if (meta::get_vol4d(parameters) *NDIM * NC * NC * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t NSPACE = parameters.get_nspace();
	int cter = 0;
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in a complex array
						hmc_complex tmp [NC][NC];
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								int pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								tmp[m][n].re = gaugefield_tmp[pos];
								tmp[m][n].im = gaugefield_tmp[pos + 1];
								cter++;
							}
						}
						//store su3matrix tmp in our format
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);

						//copy hmc_su3matrix to Matrixsu3 format
						Matrixsu3 destElem;
						destElem.e00 = tmp[0][0];
						destElem.e01 = tmp[0][1];
						destElem.e02 = tmp[0][2];
						destElem.e10 = tmp[1][0];
						destElem.e11 = tmp[1][1];
						destElem.e12 = tmp[1][2];
						destElem.e20 = tmp[2][0];
						destElem.e21 = tmp[2][1];
						destElem.e22 = tmp[2][2];

						gaugefield[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)] = destElem;
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}
}

static void copy_gaugefield_to_ildg_format(hmc_float * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters)
{
	const size_t NSPACE = parameters.get_nspace();
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);
						hmc_complex destElem [NC][NC];

						Matrixsu3 srcElem = source_in[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)];
						destElem[0][0] = srcElem.e00;
						destElem[0][1] = srcElem.e01;
						destElem[0][2] = srcElem.e02;
						destElem[1][0] = srcElem.e10;
						destElem[1][1] = srcElem.e11;
						destElem[1][2] = srcElem.e12;
						destElem[2][0] = srcElem.e20;
						destElem[2][1] = srcElem.e21;
						destElem[2][2] = srcElem.e22;

						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								dest[pos]     = destElem[m][n].re;
								dest[pos + 1] = destElem[m][n].im;
							}
						}
					}
				}
			}
		}
	}
}

hmc_float physics::lattices::Gaugefield::plaquette()
{
	assert(buffers.size() == 1);

	auto gf_dev = buffers[0];
	auto device = gf_dev->get_device();

	const hardware::buffers::Plain<hmc_float> plaq_dev(1, device);
	const hardware::buffers::Plain<hmc_float> tplaq_dev(1, device);
	const hardware::buffers::Plain<hmc_float> splaq_dev(1, device);
	device->get_gaugefield_code()->plaquette_device(gf_dev, &plaq_dev, &tplaq_dev, &splaq_dev);

	hmc_float plaq_host;
	plaq_dev.dump(&plaq_host);
	device->synchronize();
	return plaq_host;
}

static void check_sourcefileparameters(const meta::Inputparameters& parameters, const hmc_float plaquette, sourcefileparameters& parameters_source)
{
	logger.info() << "Checking sourcefile parameters against inputparameters...";
	//checking major parameters
	int tmp1, tmp2;
	std::string testobj;
	std::string msg = "Major parameters do not match: ";

	testobj = msg + "NT";
	tmp1 = parameters.get_ntime();
	tmp2 = parameters_source.lt_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NX";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.lx_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NY";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.ly_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NZ";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.lz_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "PRECISION";
	tmp1 = parameters.get_precision();
	tmp2 = parameters_source.prec_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}

	//checking minor parameters
	msg = "Minor parameters do not match: ";
	hmc_float float1, float2;
	testobj = msg + "plaquette";
	float1 = plaquette;
	float2 = parameters_source.plaquettevalue_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "beta";
	float1 = parameters.get_beta();
	float2 = parameters_source.beta_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "kappa";
	float1 = parameters.get_kappa();
	float2 = parameters_source.kappa_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "mu";
	float1 = parameters.get_mu();
	float2 = parameters_source.mu_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}

	logger.info() << "...done";
	return;
}
