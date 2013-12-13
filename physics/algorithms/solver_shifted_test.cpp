/** @file
 * Tests of the multi-shifted inverter algorithm
 * 
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "solver_shifted.hpp"
#include "rational_approximation.hpp"
#include "../../logger.hpp"
#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../lattices/staggeredfield_eo.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::solvers::solver_shifted
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(cgm_1)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//This are some possible values of sigma
	hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145,
	                                        0.12517586269872370, 1.0029328743375700};
	std::vector<hmc_float> sigma(pol, pol + sizeof(pol)/sizeof(hmc_float));
	physics::fermionmatrix::MdagM_eo matrix(system, 0.567);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	//This field is NOT that of the test explicit_stagg (D_KS_eo) because here the lattice is 4^4
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 13);
  
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-23);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	//Once obtained the solution we apply each operator (matrix + sigma) onto the field
	//out and we should obtain the field b. Hence we compare the squarenorms of b and of 
	//(matrix + sigma) * out.
	std::vector<Staggeredfield_eo*> aux;
	std::vector<hmc_float> sqnorm_out;
	hmc_float sqnorm_b = squarenorm(b);
	logger.info() << "                           sqnorm(b)=" << std::setprecision(16) << sqnorm_b;
	for(uint i=0; i<sigma.size(); i++){
		aux.push_back(new Staggeredfield_eo(system));
		matrix(aux[i], gf, *out[i]);
		saxpy(aux[i], {sigma[i],0.}, *out[i], *aux[i]);
		sqnorm_out.push_back(squarenorm(*aux[i]));
		logger.info() << "sqnorm((matrix + sigma[" << i << "]) * out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorm_b, sqnorm_out[i], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(cgm_2)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//These are some possible values of sigma
	Rational_Approximation approx(8, 1,2, 1.e-5,1);
	std::vector<hmc_float> sigma = approx.Get_b();
	physics::fermionmatrix::MdagM_eo matrix(system, 1.01335);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 123);
	
	//These are the sqnorms of the output of the CG-M algorithm from the reference code
	std::vector<hmc_float> sqnorms_ref;
	sqnorms_ref.push_back(59.877118728264179026);
	sqnorms_ref.push_back(59.875307563978779513);
	sqnorms_ref.push_back(59.865246895125224569);
	sqnorms_ref.push_back(59.811193256571264953);
	sqnorms_ref.push_back(59.521693061471694364);
	sqnorms_ref.push_back(57.988143162101970063);
	sqnorms_ref.push_back(50.327004662008008040);
	sqnorms_ref.push_back(22.236536652925686042);
	//Now I calculate the fields out
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-12);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	std::vector<hmc_float> sqnorm_out;
	for(uint i=0; i<sigma.size(); i++){
		sqnorm_out.push_back(squarenorm(*out[i]));
		logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorms_ref[i], sqnorm_out[i], 1.e-8);
	}
}

BOOST_AUTO_TEST_CASE(cgm_3)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//These are some possible values of sigma
	Rational_Approximation approx(16, 1,2, 1.e-5,1);
	
	std::vector<hmc_float> sigma = approx.Get_b();
	physics::fermionmatrix::MdagM_eo matrix(system, 0.01);
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Staggeredfield_eo b(system);
	std::vector<Staggeredfield_eo*> out;
	for(uint i=0; i<sigma.size(); i++)
		out.push_back(new Staggeredfield_eo(system));
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 123);
	
	//These are the sqnorms of the output of the CG-M algorithm from the reference code
	std::vector<hmc_float> sqnorms_ref;
	sqnorms_ref.push_back(3790.2414703421331978);
	sqnorms_ref.push_back(3789.5266994555458950);
	sqnorms_ref.push_back(3787.5285290870838253);
	sqnorms_ref.push_back(3782.6617780789301833);
	sqnorms_ref.push_back(3771.0928645597095965);
	sqnorms_ref.push_back(3743.8629449816535271);
	sqnorms_ref.push_back(3680.7256977108900173);
	sqnorms_ref.push_back(3539.0807487819220114);
	sqnorms_ref.push_back(3243.3813294554852291);
	sqnorms_ref.push_back(2708.4856234685380514);
	sqnorms_ref.push_back(1947.5005180873681638);
	sqnorms_ref.push_back(1157.0234695785266013);
	sqnorms_ref.push_back(559.56186864729193076);
	sqnorms_ref.push_back(215.02535309652577666);
	sqnorms_ref.push_back(58.529059535207736076);
	sqnorms_ref.push_back(6.2407847688851161294);
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-24);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	std::vector<hmc_float> sqnorm_out;
	for(uint i=0; i<sigma.size(); i++){
		sqnorm_out.push_back(squarenorm(*out[i]));
		logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorms_ref[i], sqnorm_out[i], 1.e-8);
	}
}

//This are just to play with cg_m to optimize it
BOOST_AUTO_TEST_CASE(cgm_4)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--nspace=4", "--ntime=4", "--fermact=rooted_stagg", "--cgmax=100000", "--cg_minimum_iteration_count=0" };
	meta::Inputparameters params(6, _params);
	hardware::System system(params);
	physics::PRNG prng(system);
	
	//These are some possible values of sigma
	Rational_Approximation approx(1, 1,2, 1.e-5,1);
	
	std::vector<hmc_float> sigma = approx.Get_b();
	physics::fermionmatrix::MdagM_eo matrix(system, 0.01);
	
	Gaugefield gf(system, prng, "hot");
	Staggeredfield_eo b(system);
 	std::vector<Staggeredfield_eo*> out;
 	for(uint i=0; i<sigma.size(); i++)
 		out.push_back(new Staggeredfield_eo(system));
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 123);
	
	int iter = cg_m(out, sigma, matrix, gf, b, system, 1.e-24);
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	std::vector<hmc_float> sqnorm_out;
	for(uint i=0; i<sigma.size(); i++){
		sqnorm_out.push_back(squarenorm(*out[i]));
		logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
	}
}

/*
BOOST_AUTO_TEST_CASE(cgm_5)
{
	using physics::fermionmatrix::DKS_eo;
	
	const char * _params[] = {"foo", "--nspace=16", "--ntime=4", "--fermact=rooted_stagg"};
	meta::Inputparameters parameters(4, _params);
	switchLogLevel(parameters.get_log_level());
	
	hardware::System system(parameters);
	physics::PRNG prng(system);
// 	physics::lattices::Gaugefield gf(system, prng);
	physics::lattices::Gaugefield gf(system, prng, "hot");
	physics::lattices::Staggeredfield_eo sf1(system);
	physics::lattices::Staggeredfield_eo sf2(system);
  
	// update gaugefield buffers once to have update links fully initialized
	gf.update_halo();
	logger.info() << "Gaugeobservables:";
	print_gaugeobservables(gf, 0);

	DKS_eo(&sf2, gf, sf1, EVEN);
	DKS_eo(&sf1, gf, sf2, ODD);
	
	for(auto dev: system.get_devices()) {
		dev->synchronize();
	}

	int hmc_iter = 2000;
	logger.info() << "Perform DKS_eo (EVEN + ODD) " << hmc_iter << " times.";
	klepsydra::Monotonic timer;
	for(int iter = 0; iter < hmc_iter; ++iter) {
		DKS_eo(&sf2, gf, sf1, EVEN);
		DKS_eo(&sf1, gf, sf2, ODD);
	}
	for(auto dev: system.get_devices()) {
		dev->synchronize();
	}
	auto elapsed_mus = timer.getTime();
		
	auto fermion_code = system.get_devices()[0]->get_fermion_staggered_code();
	size_t flop_count = fermion_code->get_flop_size("D_KS_eo");
	size_t byte_count = fermion_code->get_read_write_size("D_KS_eo");
	double gflops = static_cast<double>(flop_count) * 2 * hmc_iter / elapsed_mus / 1e3;
	double gbytes = static_cast<double>(byte_count) * 2 * hmc_iter / elapsed_mus / 1e3;
	logger.info() << "D_KS performance: " << gflops << " GFLOPS";
	logger.info() << "D_KS memory: " << gbytes << " GB/S";
	logger.info() << "Measured TIME: " << elapsed_mus / 1.e3 << "msec";
	
}
*/
