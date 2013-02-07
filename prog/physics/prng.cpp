/** @file
 * Ranlux PRNG implementation
 */

#include "prng.hpp"

#include "../host_random.h"
#include "../hardware/buffers/prng_buffer.hpp"
#include <fstream>
#include <stdexcept>
#include "../hardware/code/prng.hpp"

physics::PRNG::~PRNG()
{
for(const hardware::buffers::PRNGBuffer * buffer : buffers) {
		delete buffer;
	}
}

physics::PRNG::PRNG(const hardware::System& system) :
	system(system)
{
	using hardware::buffers::PRNGBuffer;

	auto params = system.get_inputparameters();

	// initialize host prng
	uint32_t seed = params.get_host_seed();
	prng_init(seed);

	// initialize devices
for(hardware::Device * device : system.get_devices()) {
		// create a buffer for each device
		const PRNGBuffer * buffer = new PRNGBuffer(device);
		auto code = device->get_prng_code();
		code->initialize(buffer, ++seed);
		buffers.push_back(buffer);
	}

	// additional initalization in case of known start
	if(!params.get_initial_prng_state().empty()) {
		std::ifstream file(params.get_initial_prng_state().c_str(), std::ios_base::binary);

		std::string test;
		getline(file, test);
		if(test != "OpTiMaL PRNG State") {
			throw std::invalid_argument(params.get_initial_prng_state() + " does not seem to contain a valid prng state");
		}
		file.seekg(6, std::ios_base::cur);
		size_t host_state_size = prng_size();
		int* host_state = new int[host_state_size];
		file.read(reinterpret_cast<char*>(host_state), host_state_size * sizeof(int));
		prng_set(host_state);
		delete[] host_state;
		file.seekg(1, std::ios_base::cur); // skip newline
		for(auto buffer: buffers) {
			size_t buffer_bytes;
			file >> buffer_bytes;
			if(buffer_bytes != buffer->get_bytes()) {
				throw std::invalid_argument(params.get_initial_prng_state() + " does not seem to contain a valid prng state");
			}
			file.seekg(1, std::ios_base::cur); // skip space
			char* state = new char[buffer_bytes];
			file.read(state, buffer_bytes);
			buffer->load(reinterpret_cast<const hardware::buffers::PRNGBuffer::prng_state_t *>(state));
			file.seekg(1, std::ios_base::cur); // skip newline
			delete[] state;
		}
		// TODO check if file is empty
	}
}
const std::vector<const hardware::buffers::PRNGBuffer*> physics::PRNG::get_buffers() const noexcept
{
	return buffers;
}

double physics::PRNG::get_double() noexcept {
	return prng_double();
}

Matrixsu3 physics::random_matrixsu3(physics::PRNG& prng)
{
	//this is like project_su3 from operations_gaugefield.cl, only that you do it with gaussian su3 vectors
	//taken also from the tmlqcd code
	Matrixsu3 out;
	hmc_complex a[NC];
	hmc_complex b[NC];
	hmc_complex c[NC];

	//get gaussian vector
	gaussianComplexVector(a, NC, 1., prng);
	//get norm 1
	hmc_float norm = 0;
	for(int i = 0; i < NC; i++) {
		norm += a[i].re * a[i].re + a[i].im * a[i].im;
	}
	for(int i = 0; i < NC; i++) {
		a[i].re /= sqrt(norm);
		a[i].im /= sqrt(norm);
	}

	for(;;) {
		//get another gaussian vector
		gaussianComplexVector(b, NC, 1., prng);
		//get norm 1
		norm = 0;
		for(int i = 0; i < NC; i++) {
			norm += b[i].re * b[i].re + b[i].im * b[i].im;
		}
		for(int i = 0; i < NC; i++) {
			b[i].re /= sqrt(norm);
			b[i].im /= sqrt(norm);
		}

		//calculate a* times b
		hmc_complex tmp = {0., 0.};
		for(int i = 0; i < NC; i++) {
			tmp.re += a[i].re * b[i].re + a[i].im * b[i].im;
			tmp.im += - a[i].im * b[i].re + a[i].re * b[i].im;
		}

		//project: b = tmp * a
		b[0].re -= (tmp.re * a[0].re - tmp.im * a[0].im);
		b[0].im -= (tmp.re * a[0].im + tmp.im * a[0].re);
		b[1].re -= (tmp.re * a[1].re - tmp.im * a[1].im);
		b[1].im -= (tmp.re * a[1].im + tmp.im * a[1].re);
		b[2].re -= (tmp.re * a[2].re - tmp.im * a[2].im);
		b[2].im -= (tmp.re * a[2].im + tmp.im * a[2].re);

		//get norm
		norm = 0;
		for(int i = 0; i < NC; i++) {
			norm += b[i].re * b[i].re + b[i].im * b[i].im;
		}
		norm = sqrt(norm);

		//check if norm is zero
		if (1.0 != (1.0 + norm))
			break;
	}

	// b = b/norm
	hmc_float fact = 1.0 / norm;
	b[0].re = fact * b[0].re;
	b[0].im = fact * b[0].im;
	b[1].re = fact * b[1].re;
	b[1].im = fact * b[1].im;
	b[2].re = fact * b[2].re;
	b[2].im = fact * b[2].im;

	//3rd vector, this is reconstruct12!!
	c[0].re = (a[1].re * b[2].re - a[1].im * b[2].im)
	          - (a[2].re * b[1].re - a[2].im * b[1].im);
	c[0].im = -(a[1].re * b[2].im + a[1].im * b[2].re)
	          + (a[2].re * b[1].im + a[2].im * b[1].re);

	c[1].re = (a[2].re * b[0].re - a[2].im * b[0].im)
	          - (a[0].re * b[2].re - a[0].im * b[2].im);
	c[1].im = -(a[2].re * b[0].im + a[2].im * b[0].re)
	          + (a[0].re * b[2].im + a[0].im * b[2].re);

	c[2].re = (a[0].re * b[1].re - a[0].im * b[1].im)
	          - (a[1].re * b[0].re - a[1].im * b[0].im);
	c[2].im = -(a[0].re * b[1].im + a[0].im * b[1].re)
	          + (a[1].re * b[0].im + a[1].im * b[0].re);

	//set values of matrix
	out.e00.re = a[0].re;
	out.e00.im = a[0].im;
	out.e01.re = a[1].re;
	out.e01.im = a[1].im;
	out.e02.re = a[2].re;
	out.e02.im = a[2].im;

	out.e10.re = b[0].re;
	out.e10.im = b[0].im;
	out.e11.re = b[1].re;
	out.e11.im = b[1].im;
	out.e12.re = b[2].re;
	out.e12.im = b[2].im;

	out.e20.re = c[0].re;
	out.e20.im = c[0].im;
	out.e21.re = c[1].re;
	out.e21.im = c[1].im;
	out.e22.re = c[2].re;
	out.e22.im = c[2].im;

	return out;
}

void physics::gaussianNormalPair(hmc_float * z1, hmc_float * z2, physics::PRNG& prng)
{
	hmc_float u1 = 1.0 - prng.get_double();
	hmc_float u2 = 1.0 - prng.get_double();
	hmc_float p  = sqrt(-2 * log(u1));
	*z1 = p * cos(2 * PI * u2);
	*z2 = p * sin(2 * PI * u2);
	// SL: not yet tested
}

void physics::gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma, physics::PRNG& prng)
{
	// SL: this fills real and imaginary part of a vector of "length" complex numbers
	//     with components drawn with a Gaussian distribution and variance sigma
	for(int idx = 0; idx < length; idx++) {
		gaussianNormalPair(&vector[idx].re, &vector[idx].im, prng);
		vector[idx].re *= sigma;
		vector[idx].im *= sigma;
	}
	// SL: not yet tested
}

void physics::PRNG::store(const std::string filename) const
{
	// TODO this misses a lot of error handling
	std::ofstream file(filename.c_str(), std::ios_base::binary);
	file << "OpTiMaL PRNG State\n";
	file << "Host: ";
	size_t host_state_size = prng_size();
	int* host_state = new int[host_state_size];
	prng_get(host_state);
	file.write(reinterpret_cast<char*>(host_state), host_state_size * sizeof(int));
	delete[] host_state;
	file << '\n';
	for(auto buffer: buffers) {
		size_t buffer_bytes = buffer->get_bytes();
		file << buffer_bytes << ' ';
		char* state = new char[buffer_bytes];
		buffer->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(state));
		file.write(state, buffer_bytes);
		file << '\n';
		delete[] state;
	}
}

