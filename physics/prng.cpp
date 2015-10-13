/** @file
 * Ranlux PRNG implementation
 *
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "prng.hpp"

#include "../host_functionality/host_random.h"
#include "../hardware/buffers/prng_buffer.hpp"
#include <fstream>
#include <stdexcept>
#include "../hardware/device.hpp"
#include "../hardware/code/prng.hpp"
#include "utilities.hpp"

//http://www.ridgesolutions.ie/index.php/2013/05/30/boost-link-error-undefined-reference-to-boostfilesystemdetailcopy_file/
#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>

physics::PRNG::~PRNG()
{
	for(const hardware::buffers::PRNGBuffer * buffer : buffers)
	{
		delete buffer;
	}
}

void readFirstLine_checkForTag( std::ifstream & file, const std::string filename)
{
	std::string test;
	const std::string tag = "OpTiMaL PRNG State";
	getline(file, test);
	if(test != tag) {
		logger.fatal() << "Did not find correct tag in prng-file, found \""<< test << "\" instead!";
		throw std::invalid_argument("\"" + filename + "\" does not seem to contain a valid prng state");
	}
}

int* readSecondLine_extractHostSeed( std::ifstream & file )
{
	file.seekg(6, std::ios_base::cur);
	size_t host_state_size = prng_size();
	int* host_state = new int[host_state_size];
	file.read(reinterpret_cast<char*>(host_state), host_state_size * sizeof(int));
	file.seekg(1, std::ios_base::cur); // skip newline
	return host_state;
}

struct PrngStateContent
{
	PrngStateContent( size_t bufferBytesIn, std::unique_ptr<char[]> contentIn) :
		bufferBytes( bufferBytesIn ), content(std::move(contentIn)) {}
	size_t bufferBytes;
	std::unique_ptr<char[]> content;
};

PrngStateContent readLine_prngState( std::ifstream & file )
{
	size_t buffer_bytes;
	file >> buffer_bytes;
	file.seekg(1, std::ios_base::cur); // skip space
    std::unique_ptr<char[]> state(new char[buffer_bytes]);
	file.read(state.get(), buffer_bytes);
	file.seekg(1, std::ios_base::cur); // skip newline

	return PrngStateContent(buffer_bytes, std::move(state));
}

struct PrngFileReader
{
	PrngFileReader( std::string filename, const uint numberOfExpectedBuffers )
	{
		std::ifstream file(filename.c_str(), std::ios_base::binary);
		readFirstLine_checkForTag( file, filename );
		hostState = readSecondLine_extractHostSeed( file );
		for( uint currentBuffer = 0; currentBuffer < numberOfExpectedBuffers; currentBuffer ++)
		{
			prngStates.push_back( readLine_prngState( file ) );
		}
	}
	~PrngFileReader()
	{
		delete[] hostState;
	}
	std::vector<PrngStateContent> prngStates;
	int * hostState;
};

physics::PRNG::PRNG(const hardware::System& system, const physics::ParametersPrngInterface * parametersIn) :
	system(system), parameters(parametersIn)
{
	using hardware::buffers::PRNGBuffer;

	// initialize host prng
	uint32_t seed = parameters->getHostSeed();
	prng_init(seed);

	// initialize devices
	for(hardware::Device * device : system.get_devices()) {
		// create a buffer for each device
		const PRNGBuffer * buffer = new PRNGBuffer(device, parameters->useSameRandomNumbers());
		auto code = device->getPrngCode();
		code->initialize(buffer, ++seed);
		buffers.push_back(buffer);
	}

	// additional initialization in case of known start
	if(!parameters->getInitialPrngStateFilename().empty())
	{
		logger.debug() << "Read prng state from file \"" + parameters->getInitialPrngStateFilename() + "\"...";
		PrngFileReader fileContent( parameters->getInitialPrngStateFilename(), buffers.size() );

		prng_set(fileContent.hostState);
		for( uint currentBufferIndex = 0; currentBufferIndex < buffers.size(); currentBufferIndex ++ )
		{
			size_t buffer_bytes = fileContent.prngStates.at(currentBufferIndex).bufferBytes;
			if(buffer_bytes != buffers.at(currentBufferIndex)->get_bytes()) {
				throw std::invalid_argument(parameters->getInitialPrngStateFilename() + " does not seem to contain a valid prng state");
			}
			buffers.at(currentBufferIndex)->load(reinterpret_cast<const hardware::buffers::PRNGBuffer::prng_state_t *>( fileContent.prngStates.at(currentBufferIndex).content.get() ));
		}
	}
	else
	{
		logger.debug() << "Did not find prng file...";
	}
}
const std::vector<const hardware::buffers::PRNGBuffer*> physics::PRNG::get_buffers() const noexcept
{
	return buffers;
}

double physics::PRNG::get_double() const noexcept{
	return prng_double();
}

Matrixsu3 physics::random_matrixsu3(const physics::PRNG& prng)
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

void physics::gaussianNormalPair(hmc_float * z1, hmc_float * z2, const physics::PRNG& prng)
{
	hmc_float u1 = 1.0 - prng.get_double();
	hmc_float u2 = 1.0 - prng.get_double();
	hmc_float p  = sqrt(-2 * log(u1));
	*z1 = p * cos(2 * PI * u2);
	*z2 = p * sin(2 * PI * u2);
	// SL: not yet tested
}

void physics::gaussianComplexVector(hmc_complex * vector, int length, hmc_float sigma, const physics::PRNG& prng)
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

void physics::PRNG::verifyWritingWasSuccessful(const std::string filename) const
{
	logger.info() << "Verify writing to file: Re-read prng state from file \"" + filename + "\"...";
	PrngFileReader fileContent( filename, buffers.size() );

	for(size_t i = 0; i < this->get_buffers().size(); ++i) {
		auto buf1 = this->get_buffers().at(i);

		char* prng1_state = new char[buf1->get_bytes()];
		buf1->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng1_state));

		auto prng2_state = fileContent.prngStates.at(i).content.get();

		for( size_t j = 0; j < buf1->get_bytes(); j++)
		{
			//logger.fatal() << prng1_state[j] << " " << prng2_state[j] << " " << (prng1_state[j] == prng2_state[j]);
			if (prng1_state[j] != prng2_state[j])
			{
				logger.fatal() << "Writing of prng file not successful! Aborting...";
				throw File_Exception(filename);
			}
		}

		delete[] prng1_state;
	}
	logger.info() << "...done";
}

void physics::PRNG::store(const std::string filename) const
{
	logger.info() << "saving current prng state to file \"" << filename << "\"";
	if ( boost::filesystem::exists( filename ) )
	{
	  //TODO: It seems that the function boost::filesystem::copy_file gives a linking error in compilation
	  //      /cm/shared/apps/boost/1.52.0/include/boost/filesystem/operations.hpp:381: undefined reference to `boost::filesystem::detail::copy_file(...)
	  //      therefore here we rename files and then delete. This should nevertheless be fixed!
	  const std::string backupFilename = filename + "_backup";
	  if(boost::filesystem::exists(backupFilename))
	    boost::filesystem::rename(backupFilename, backupFilename + "_tmp");
	  logger.warn() << "Found existing file of name \"" << filename << "\". Store this to \"" << backupFilename << "\"...";
	  boost::filesystem::rename( filename, backupFilename);
	  if(boost::filesystem::exists(backupFilename + "_tmp"))
	    boost::filesystem::remove(backupFilename + "_tmp");
	  //boost::filesystem::copy_file( filename, backupFilename, boost::filesystem::copy_option::overwrite_if_exists);
	}

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
	file.close();

	verifyWritingWasSuccessful( filename );
}

std::string physics::PRNG::getName(int number) const noexcept
{
	return physics::getConfigurationName( parameters->getNamePrefix(), parameters->getNamePostfix(), parameters->getNumberOfDigitsInName(), number);
}

void physics::PRNG::save()
{
	store(this->getName());
}

void physics::PRNG::saveToSpecificFile(int number)
{
	store(this->getName(number));
}

bool physics::PRNG::operator == (const physics::PRNG & prng) const
{
	for(size_t i = 0; i < this->get_buffers().size(); ++i) {
		auto buf1 = this->get_buffers().at(i);
		auto buf2 = prng.get_buffers().at(i);

		char* prng1_state = new char[buf1->get_bytes()];
		buf1->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng1_state));

		char* prng2_state = new char[buf2->get_bytes()];
		buf2->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng2_state));

		for( size_t j = 0; j < buf1->get_bytes(); j++)
		{
			if (prng1_state[j] != prng2_state[j])
			{
				return false;
			}
		}

		delete[] prng1_state;
		delete[] prng2_state;
	}
	return true;
}

bool physics::PRNG::operator != (const physics::PRNG & prng) const
{
	return *this == prng ? false : true;
}
