/** @file
 *
 * Copyright (c) 2014,2018 Alessandro Sciarra
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
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

#include "parametersObs.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

static common::pbp_version translatePbpVersionToEnum(std::string);

bool meta::ParametersObs::get_measure_transportcoefficient_kappa() const noexcept
{
    return measure_transportcoefficient_kappa;
}

bool meta::ParametersObs::get_measure_rectangles() const noexcept
{
    return measure_rectangles;
}

bool meta::ParametersObs::get_measure_correlators() const noexcept
{
    return measure_correlators;
}

bool meta::ParametersObs::get_measure_pbp() const noexcept
{
    return measure_pbp;
}

common::pbp_version meta::ParametersObs::get_pbp_version() const noexcept
{
    return pbp_version_;
}

int meta::ParametersObs::get_corr_dir() const noexcept
{
    return corr_dir;
}

int meta::ParametersObs::get_pbp_measurements() const noexcept
{
    return pbp_measurements;
}

meta::ParametersObs::ParametersObs()
    : measure_transportcoefficient_kappa(false)
    , measure_rectangles(false)
    , measure_correlators(false)
    , measure_pbp(false)
    , corr_dir(3)
    , pbp_measurements(1)
    , options("Observables options")
    , pbp_version_String("std")
    , pbp_version_(common::pbp_version::std)
{
    // clang-format off
    options.add_options()
    ("correlatorDirection", po::value<int>(&corr_dir)->default_value(corr_dir), "The direction for the correlator.")
    ("measureCorrelators", po::value<bool>(&measure_correlators)->default_value(measure_correlators), "Whether to measure fermionic correlators.")
    ("measurePbp", po::value<bool>(&measure_pbp)->default_value(measure_pbp), "Whether to measure chiral condensate.")
    ("pbpVersion",  po::value<std::string>(&pbp_version_String)->default_value(pbp_version_String), "Which version of chiral condensate to measure (one among 'std' and 'tm_one_end_trick').")
    ("pbpMeasurements", po::value<int>(&pbp_measurements)->default_value(pbp_measurements), "Number of chiral condensate measurements (for 'rooted_stagg' fermion action only!).")
    ("measureTransportCoefficientKappa", po::value<bool>(&measure_transportcoefficient_kappa)->default_value(measure_transportcoefficient_kappa), "Whether to measure the transport coefficient kappa.")
    ("measureRectangles", po::value<bool>(&measure_rectangles)->default_value(measure_rectangles), "Whether to measure rectangles.");
    // clang-format on
}

static common::pbp_version translatePbpVersionToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::pbp_version> m;
    m["std"]              = common::std;
    m["tm_one_end_trick"] = common::tm_one_end_trick;

    common::pbp_version a = m[s];
    if (a) {  // map returns 0 if element is not found
        return a;
    } else {
        throw Invalid_Parameters("Invalid pbp version!", "std, tm_one_end_trick", s);
    }
}

void meta::ParametersObs::makeNeededTranslations()
{
    pbp_version_ = translatePbpVersionToEnum(pbp_version_String);
}
