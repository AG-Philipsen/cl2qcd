/** @file
 * Tests of the flavour doublets algorithms
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
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

#include "staggeredChiralCondensate.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::chiral_condensate_stagg
#include <boost/test/unit_test.hpp>

#include "../prng.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

/* Here pbp_ref_im_minmax are the minimum and maximum pbp immaginary part obtained in the reference
 * code in 100 measurements. This is done because of big fluctuations: a check of the closeness of
 * the values up to 300% or even more is nonsense.
 *  - pbp_ref_re is the mean of the ref code 100 measurements
 *  - pbp_ref_im_minmax.im is maximum
 *  - pbp_ref_im_minmax.re is minimum
 */
static void test_chiral_condensate_stagg(std::string content, hmc_float pbp_ref_re, hmc_complex pbp_ref_im_minmax, bool cold, bool thetaT)
{
	using namespace physics::lattices;

	std::vector<const char*> options(1, "foo");
	options.push_back("--ntime=4");
	options.push_back("--num_dev=1");
	options.push_back("--fermact=rooted_stagg");
	options.push_back("--num_tastes=2");
	options.push_back("--mass=0.1");
	options.push_back("--theta_fermion_spatial=0");
	if(thetaT)
	  options.push_back("--theta_fermion_temporal=1");
	else
	  options.push_back("--theta_fermion_temporal=0");
	options.push_back("--num_sources=10");
 	options.push_back("--sourcetype=volume");
	std::string tmp = "--sourcecontent=" + content;
	options.push_back(tmp.c_str());
	
	meta::Inputparameters params(11, &(options[0]));
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	const Gaugefield *gf;
	if(cold)
	  gf = new Gaugefield(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
	else //This configuration for the Ref.Code is the same as for example dks_input_5
	  gf = new Gaugefield(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	
	hmc_complex pbp = physics::observables::staggered::measureChiralCondensate(*gf, prng, system, interfacesHandler);
	
	logger.info() << "Chiral condensate pbp = (" << std::setprecision(12) << pbp.re << ", " << pbp.im << ")";
	if(content == "one"){
	  BOOST_CHECK_CLOSE(pbp.re, pbp_ref_re, 1.e-8);
	  //For some reason, the result could not be equal to the number 0, but below the double precision
	  //(on the Loewe this is the case). Then we check this to make the test pass (for the machine
	  //all numbers below 1.e-16 are in practice zero).
	  if(pbp_ref_im_minmax.re == 0)
	      BOOST_CHECK_SMALL(pbp.im, 1.e-15);
	  else
	      BOOST_CHECK_CLOSE(pbp.im, pbp_ref_im_minmax.re, 1.e-8);
	  if(pbp_ref_im_minmax.im == 0)
	      BOOST_CHECK_SMALL(pbp.im, 1.e-15);
	  else
	      BOOST_CHECK_CLOSE(pbp.im, pbp_ref_im_minmax.im, 1.e-8);
	}else{
	  BOOST_CHECK_CLOSE(pbp.re, pbp_ref_re, 25); //25% is needed because rand num. are not the same in ref code
	  BOOST_CHECK_PREDICATE( std::less_equal<hmc_float>(), (pbp.im)(pbp_ref_im_minmax.im) );
	  BOOST_CHECK_PREDICATE( std::greater_equal<hmc_float>(), (pbp.im)(pbp_ref_im_minmax.re) );
	}
}


BOOST_AUTO_TEST_CASE(volume_source_one)
{
	test_chiral_condensate_stagg("one", 15.0, {0.0, 0.0}, true, 0);
	test_chiral_condensate_stagg("one", 0.29411764705882464943, {0.0, 0.0}, true, 1);
	test_chiral_condensate_stagg("one", 0.20344039707296779351, 
				     {-0.0065506178486867301658, -0.0065506178486867301658}, false, 0);
	test_chiral_condensate_stagg("one", 0.25968115501530797395,
				     {-0.021963576512093983817, -0.021963576512093983817}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_z4)
{
	test_chiral_condensate_stagg("z4", 1.0194940543945434364,
				     {-0.030501775469794135953, 0.024558907459382613159}, true, 0);
 	test_chiral_condensate_stagg("z4", 0.10173717281774301291,
				     {-0.05257050359883207874, 0.038588374094059074704}, true, 1);
	test_chiral_condensate_stagg("z4", 0.21838899047431276079,
				     {-0.051483155453766353549, 0.046625129123589938163}, false, 0);
	test_chiral_condensate_stagg("z4", 0.25388166919474708383,
				     {-0.055589254436301166473, 0.059235978143224822523}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_gaussian)
{
	test_chiral_condensate_stagg("gaussian", 1.0125526769942891914,
				     {-0.037359420873174425948, 0.029119643958998853162}, true, 0);
	test_chiral_condensate_stagg("gaussian", 0.10181163880167679037,
				     {-0.033993182981299545353, 0.030919993128422099127}, true, 1);
	test_chiral_condensate_stagg("gaussian", 0.21950471423558026718,
				     {-0.058332954616782443924, 0.042412447314374657203}, false, 0);
	test_chiral_condensate_stagg("gaussian", 0.25382073751894007607,
				     {-0.06509237221358404879, 0.050458791308785347351}, false, 1);
}

BOOST_AUTO_TEST_CASE(volume_source_z2)
{
	test_chiral_condensate_stagg("z2", 1.0215954074131750051, {0.0, 0.0}, true, 0);
	test_chiral_condensate_stagg("z2", 0.10167935377929424035,
				     {-0.023941038345876981125, 0.024737496665225969239}, true, 1);
	test_chiral_condensate_stagg("z2", 0.21903384850763127356, 
				     {-0.051079063817292484628, 0.051362451897559530112}, false, 0);
	test_chiral_condensate_stagg("z2", 0.25440057886149292088,
				     {-0.045645403465476366844, 0.051219156998163921368}, false, 1);
}
