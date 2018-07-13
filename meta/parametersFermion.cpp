/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#include "parametersFermion.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

static common::action translateFermionActionToEnum(std::string);

bool meta::ParametersFermion::get_use_chem_pot_re() const noexcept
{
    return use_chem_pot_re;
}
bool meta::ParametersFermion::get_use_chem_pot_im() const noexcept
{
    return use_chem_pot_im;
}

common::action meta::ParametersFermion::get_fermact() const noexcept
{
    return fermact;
}
common::action meta::ParametersFermion::get_fermact_mp() const noexcept
{
    return fermactMP;
}
double meta::ParametersFermion::get_kappa() const noexcept
{
    return kappa;
}
double meta::ParametersFermion::get_mass() const noexcept
{
    return mass;
}
double meta::ParametersFermion::get_mu() const noexcept
{
    return mu;
}
double meta::ParametersFermion::get_csw() const noexcept
{
    return csw;
}
double meta::ParametersFermion::get_kappa_mp() const noexcept
{
    return kappa_mp;
}
double meta::ParametersFermion::get_mu_mp() const noexcept
{
    return mu_mp;
}
double meta::ParametersFermion::get_csw_mp() const noexcept
{
    return csw_mp;
}

double meta::ParametersFermion::get_theta_fermion_spatial() const noexcept
{
    return theta_fermion_spatial;
}
double meta::ParametersFermion::get_theta_fermion_temporal() const noexcept
{
    return theta_fermion_temporal;
}
double meta::ParametersFermion::get_chem_pot_re() const noexcept
{
    return chem_pot_re;
}
double meta::ParametersFermion::get_chem_pot_im() const noexcept
{
    return chem_pot_im;
}
bool meta::ParametersFermion::get_use_eo() const noexcept
{
    return use_eo;
}

bool meta::ParametersFermion::get_use_merge_kernels_fermion() const noexcept
{
    return use_merge_kernels_fermion;
}
bool meta::ParametersFermion::get_use_merge_kernels_spinor() const noexcept
{
    return use_merge_kernels_spinor;
}

meta::ParametersFermion::ParametersFermion()
    : kappa(0.125)
    , mass(0.1)
    , mu(0.006)
    , csw(0.)
    , kappa_mp(0.125)
    , mu_mp(0.006)
    , csw_mp(0.)
    , theta_fermion_spatial(0.)
    , theta_fermion_temporal(0.)
    , use_chem_pot_re(false)
    , use_chem_pot_im(false)
    , chem_pot_re(0.)
    , chem_pot_im(0.)
    , use_eo(true)
    , use_merge_kernels_fermion(false)
    , use_merge_kernels_spinor(false)
    , options("Fermion options")
    , fermactString("wilson")
    , fermactMPString("wilson")
    , fermact(common::action::wilson)
    , fermactMP(common::action::wilson)
{
    // clang-format off
    options.add_options()
    ("fermionAction", po::value<std::string>(&fermactString)->default_value(fermactString),"Which type of fermion action to use (e.g. wilson, twistedmass, rooted_stagg).")
    ("fermionActionMP", po::value<std::string>(&fermactMPString)->default_value(fermactMPString),"Which type of fermion action to use in the Mass Preconditioning trick when it is switched on (e.g. wilson).")
    //todo: change this default value!
    ("kappa", po::value<double>(&kappa)->default_value(kappa),"The hopping parameter in in the 'wilson' action.")
    ("mass", po::value<double>(&mass)->default_value(mass),"The bare quark mass in the 'rooted_stagg' action.")
    ("mu", po::value<double>(&mu)->default_value(mu),"The twisted mass parameter in the 'twistedmass' action.")
    ("csw", po::value<double>(&csw)->default_value(csw),"The clover coefficient in the 'clover' action.")
    ("kappaMP", po::value<double>(&kappa_mp)->default_value(kappa_mp),"The hopping parameter in the 'wilson' action part with Mass Preconditioning.")
    ("muMP", po::value<double>(&mu_mp)->default_value(mu_mp),"The twisted mass parameter in the 'twistedmass' action part with Mass Preconditioning.")
    ("cswMP", po::value<double>(&csw_mp)->default_value(csw_mp),"The clover coefficient in the 'clover' action part with Mass Preconditioning.")
    ("thetaFermionSpatial", po::value<double>(&theta_fermion_spatial)->default_value(theta_fermion_spatial),"The fermion boundary condition phase in spatial direction (e.g. 0 or 1 for periodic or anti-periodic BC, respectively).")
    ("thetaFermionTemporal", po::value<double>(&theta_fermion_temporal)->default_value(theta_fermion_temporal),"The fermion boundary condition phase in temporal direction (e.g. 0 or 1 for periodic or anti-periodic BC, respectively).")
    ("useChemicalPotentialRe", po::value<bool>(&use_chem_pot_re)->default_value(use_chem_pot_re),"Whether to switch on a nonzero real part of the quark chemical potential.")
    ("useChemicalPotentialIm", po::value<bool>(&use_chem_pot_im)->default_value(use_chem_pot_im),"Whether to switch on a nonzero imaginary part of the quark chemical potential.")
    ("chemicalPotentialRe", po::value<double>(&chem_pot_re)->default_value(chem_pot_re),"The value of the real part of the quark chemical potential.")
    ("chemicalPotentialIm", po::value<double>(&chem_pot_im)->default_value(chem_pot_im),"The value of the imaginary part of the quark chemical potential.")
    ("useEO", po::value<bool>(&use_eo)->default_value(use_eo),"Whether to switch on Even Odd preconditioning.")
    ("useKernelMergingSpinor", po::value<bool>(&use_merge_kernels_spinor)->default_value(use_merge_kernels_spinor), "Whether to use kernel merging for spinor kernels.")
    ("useKernelMergingFermionMatrix", po::value<bool>(&use_merge_kernels_fermion)->default_value(use_merge_kernels_fermion), "Whether to use kernel merging for fermion matrix kernels.");
    // clang-format on
}

static common::action translateFermionActionToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::action> m;
    m["wilson"]       = common::action::wilson;
    m["clover"]       = common::action::clover;
    m["twistedmass"]  = common::action::twistedmass;
    m["rooted_stagg"] = common::action::rooted_stagg;

    common::action a = m[s];
    if (a) {  // map returns 0 if element is not found
        return a;
    } else {
        throw Invalid_Parameters("Unkown fermion action!", "wilson, clover, twistedmass, rooted_stagg", s);
    }
}

void meta::ParametersFermion::makeNeededTranslations()
{
    fermact   = translateFermionActionToEnum(fermactString);
    fermactMP = translateFermionActionToEnum(fermactMPString);
}
