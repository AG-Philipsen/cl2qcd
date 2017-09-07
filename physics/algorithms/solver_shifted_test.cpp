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
#include "../../host_functionality/logger.hpp"
#include "../lattices/util.hpp"
#include "../../meta/util.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::solvers::solver_shifted
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(cgm_1)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	
	meta::Inputparameters* params;
	for(int i=0; i<2; i++){
	    if(i==0){
	        const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.567", "--num_dev=1"};
	        params = new meta::Inputparameters(5, _params);
	    }else{
	        const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.567", "--use_merge_kernels_spinor=true", "--num_dev=1"};
	        params = new meta::Inputparameters(6, _params);
	    }

	    hardware::HardwareParametersImplementation hP(params);
	    hardware::code::OpenClKernelParametersImplementation kP(*params);
	    hardware::System system(hP, kP);
	    physics::InterfacesHandlerImplementation interfacesHandler{*params};
	    physics::PrngParametersImplementation prngParameters{*params};
	    physics::PRNG prng{system, &prngParameters};

	    //This are some possible values of sigma
	    hmc_float pol[5] = {0.0002065381736724, 0.00302707751065980, 0.0200732678058145,
	            0.12517586269872370, 1.0029328743375700};
	    std::vector<hmc_float> sigma(pol, pol + sizeof(pol)/sizeof(hmc_float));
	    physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());

	    //This configuration for the Ref.Code is the same as for example dks_input_5
	    Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	    Staggeredfield_eo b(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

        std::vector<std::shared_ptr<Staggeredfield_eo> > out;
	    for(uint i=0; i<sigma.size(); i++)
	        out.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));

	    //This field is NOT that of the test explicit_stagg (D_KS_eo) because here the lattice is 4^4
	    pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 13);

	    int iter = cg_m(out, matrix, gf, sigma, b, system, interfacesHandler, 1.e-23, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	    logger.info() << "CG-M algorithm converged in " << iter << " iterations.";

	    //Once obtained the solution we apply each operator (matrix + sigma) onto the field
	    //out and we should obtain the field b. Hence we compare the squarenorms of b and of
	    //(matrix + sigma) * out.
	    std::vector<std::shared_ptr<Staggeredfield_eo> > aux;
	    std::vector<hmc_float> sqnorm_out;
	    hmc_float sqnorm_b = squarenorm(b);
	    logger.info() << "                           sqnorm(b)=" << std::setprecision(16) << sqnorm_b;
	    for(uint i=0; i<sigma.size(); i++){
	        aux.push_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
	        matrix(aux[i].get(), gf, *out[i], &(interfacesHandler.getAdditionalParameters<Staggeredfield_eo>()));
	        saxpy(aux[i].get(), {sigma[i],0.}, *out[i], *aux[i]);
	        sqnorm_out.push_back(squarenorm(*aux[i]));
	        logger.info() << "sqnorm((matrix + sigma[" << i << "]) * out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
	        BOOST_CHECK_CLOSE(sqnorm_b, sqnorm_out[i], 1.e-8);
	    }
	    delete params;
	}
}

BOOST_AUTO_TEST_CASE(cgm_2)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	meta::Inputparameters* params;
	for(int i=0; i<2; i++){
	    if(i==0){
	        const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=1.01335", "--num_dev=1"};
	        params = new meta::Inputparameters(5, _params);
	    }else{
	        const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=1.01335", "--use_merge_kernels_spinor=true", "--num_dev=1"};
	        params = new meta::Inputparameters(6, _params);
	    }

        hardware::HardwareParametersImplementation hP(params);
        hardware::code::OpenClKernelParametersImplementation kP(*params);
        hardware::System system(hP, kP);
	    physics::InterfacesHandlerImplementation interfacesHandler{*params};
	    physics::PrngParametersImplementation prngParameters{*params};
	    physics::PRNG prng{system, &prngParameters};

	    //These are some possible values of sigma
	    Rational_Approximation approx(8, 1,2, 1.e-5,1);
	    std::vector<hmc_float> sigma = approx.Get_b();
	    physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());

	    //This configuration for the Ref.Code is the same as for example dks_input_5
	    Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	    Staggeredfield_eo b(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	    std::vector<std::shared_ptr<Staggeredfield_eo> > out;
	    for(uint i=0; i<sigma.size(); i++)
	        out.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
	    //This field is that of the test explicit_stagg, part 2 (D_KS_eo) (ref_vec_odd because the seed is 123)
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
        int iter = cg_m(out, matrix, gf, sigma, b, system, interfacesHandler, 1.e-12, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	    logger.info() << "CG-M algorithm converged in " << iter << " iterations.";

	    std::vector<hmc_float> sqnorm_out;
	    for(uint i=0; i<sigma.size(); i++){
	        sqnorm_out.push_back(squarenorm(*out[i]));
	        logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
	        BOOST_CHECK_CLOSE(sqnorms_ref[i], sqnorm_out[i], 1.e-8);
	    }
	    delete params;
	}
}

BOOST_AUTO_TEST_CASE(cgm_3)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.01", "--num_dev=1"};
	meta::Inputparameters params(5, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	
	//These are some possible values of sigma
	Rational_Approximation approx(16, 1,2, 1.e-5,1);
	
	std::vector<hmc_float> sigma = approx.Get_b();
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Staggeredfield_eo b(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    std::vector<std::shared_ptr<Staggeredfield_eo> > out;
    for(uint i=0; i<sigma.size(); i++)
        out.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo) (ref_vec_odd because the seed is 123)
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
	int iter = cg_m(out, matrix, gf, sigma, b, system, interfacesHandler, 1.e-24, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	std::vector<hmc_float> sqnorm_out;
	for(uint i=0; i<sigma.size(); i++){
		sqnorm_out.push_back(squarenorm(*out[i]));
		logger.info() << "sqnorm(out[" << i << "])=" << std::setprecision(16) << sqnorm_out[i];
		BOOST_CHECK_CLOSE(sqnorms_ref[i], sqnorm_out[i], 1.e-8);
	}
}

//This test is to test if the cgm works correctly in the case it is used
//as a standard CG. This means that here we set the shift to zero and the
//number if equations to 1. In the reference code we will produce the result with
//the standard CG, because it is there available. The reason for this test is that,
//until a standard CG for staggered field will be available, the cgm will be
//used in the chiral condensate evaluation as a standard CG.
BOOST_AUTO_TEST_CASE(cgm_4)
{
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	using namespace physics::algorithms;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.01", "--num_dev=1"};
	meta::Inputparameters params(5, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};

	std::vector<hmc_float> sigma(1, 0.0);
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	Staggeredfield_eo b(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    std::vector<std::shared_ptr<Staggeredfield_eo> > out;
    for(uint i=0; i<sigma.size(); i++)
        out.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
	
	{
	//Cold b
	b.set_cold();
	//These is the sqnorm of the output of the CG algorithm from the reference code
	hmc_float sqnorms_ref = 9.0597433493689383255;
	int iter = cg_m(out, matrix, gf, sigma, b, system, interfacesHandler, 1.e-24, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	hmc_float sqnorm_out = squarenorm(*out[0]);
	logger.info() << "sqnorm(out)=" << std::setprecision(16) << sqnorm_out;
	BOOST_CHECK_CLOSE(sqnorms_ref, sqnorm_out, 1.e-8);
	}
	
	{
	//This field is that of the test explicit_stagg, part 2 (D_KS_eo) (ref_vec_odd because the seed is 123)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&b, 123);
	//These is the sqnorm of the output of the CG algorithm from the reference code
	hmc_float sqnorms_ref = 3790.3193634090343949;
	int iter = cg_m(out, matrix, gf, sigma, b, system, interfacesHandler, 1.e-24, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	logger.info() << "CG-M algorithm converged in " << iter << " iterations.";
	
	hmc_float sqnorm_out = squarenorm(*out[0]);
	logger.info() << "sqnorm(out)=" << std::setprecision(16) << sqnorm_out;
	BOOST_CHECK_CLOSE(sqnorms_ref, sqnorm_out, 1.e-8);
	}
}

