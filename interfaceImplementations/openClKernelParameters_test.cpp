/**
 * Copyright 2015 Francesca Cuteri, christopher Pinke
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

#include "openClKernelParameters.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::code::openClKernelParameters
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(implementByMeansOfMetaInputparameters)
{
	const char * argv []  = {"foo"};
	const meta::Inputparameters fullParameters{1, argv};
	hardware::code::OpenClKernelParametersImplementation openClKernelParameters( fullParameters );

	BOOST_REQUIRE_EQUAL( openClKernelParameters.getNs() , fullParameters.get_nspace() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getNt() , fullParameters.get_ntime() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getPrecision() , fullParameters.get_precision() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseChemPotRe(), fullParameters.get_use_chem_pot_re() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getChemPotRe(), fullParameters.get_chem_pot_re() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseChemPotIm(), fullParameters.get_use_chem_pot_im() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getChemPotIm(), fullParameters.get_chem_pot_im() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseSmearing(), fullParameters.get_use_smearing() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getRho(), fullParameters.get_rho() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getRhoIter(), fullParameters.get_rho_iter() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseRec12(), fullParameters.get_use_rec12() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseEo(), fullParameters.get_use_eo() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getFermact(), fullParameters.get_fermact() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMetroApproxOrd(), fullParameters.get_metro_approx_ord() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMdApproxOrd(), fullParameters.get_md_approx_ord() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getThetaFermionSpatial(), fullParameters.get_theta_fermion_spatial() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getThetaFermionTemporal(), fullParameters.get_theta_fermion_temporal() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getBeta(), fullParameters.get_beta() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getKappa(), fullParameters.get_kappa() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getNumSources(), fullParameters.get_num_sources() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getSourceContent(), fullParameters.get_sourcecontent() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseAniso(), fullParameters.get_use_aniso() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getSpatialLatticeVolume(), meta::get_volspace(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getLatticeVolume(), meta::get_vol4d(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseRectangles(), meta::get_use_rectangles(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getC0(), meta::get_c0(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getC1(), meta::get_c1(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getXi0(), meta::get_xi_0(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getFloatSize(), meta::get_float_size(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMatSize(), meta::get_mat_size(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getSpinorFieldSize(), meta::get_vol4d(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getEoprecSpinorFieldSize(), meta::get_vol4d(fullParameters)/2 );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getCorrDir(), fullParameters.get_corr_dir() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMeasureCorrelators(), fullParameters.get_measure_correlators() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseMergeKernelsFermion(), fullParameters.get_use_merge_kernels_fermion() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseMergeKernelsSpinor(), fullParameters.get_use_merge_kernels_spinor() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMuBar(), meta::get_mubar(fullParameters) );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getMass(), fullParameters.get_mass() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getUseSameRndNumbers(), fullParameters.get_use_same_rnd_numbers() );
	BOOST_REQUIRE_EQUAL( openClKernelParameters.getHostSeed(), fullParameters.get_host_seed() );
}
