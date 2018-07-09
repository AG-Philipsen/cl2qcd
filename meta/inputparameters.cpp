/** @file
 * Input file handling implementation
 *
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2012,2014,2015 Christopher Pinke
 * Copyright (c) 2014,2015,2018 Alessandro Sciarra
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

#include "inputparameters.hpp"

#include "config_file_normalizer.hpp"

#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using namespace meta;
/**
 * Adds all alternative option names to the ConfigFileNormlizer instance
 */
static void add_option_aliases(meta::ConfigFileNormalizer* const);

void meta::Inputparameters::ChecksStringOptionsAndMapToEnum()
{
    ParametersGauge::makeNeededTranslations();
    ParametersConfig::makeNeededTranslations();
    ParametersFermion::makeNeededTranslations();
    ParametersIntegrator::makeNeededTranslations();
    ParametersObs::makeNeededTranslations();
    ParametersSolver::makeNeededTranslations();
    ParametersSources::makeNeededTranslations();
}

Inputparameters::Inputparameters(int argc, const char** argv, std::string parameterSet)
{
    logger.debug() << "Read in parameters...";
    /**
     * First handle all the stuff that can only be done on the cmd-line.
     * We need that to get the option file.
     */
    meta::InputparametersOptions cmd_opts("Generic options");
    // clang-format off
    cmd_opts.add_options()
    ("help,h", "Produce helper for the specific executable.")
    ("inputFile", po::value<std::string>(), "The path of the file containing the input parameters.");
    // clang-format on
    // TODO add log-level etc
    po::positional_options_description pos_opts;
    pos_opts.add("inputFile", 1);

    meta::InputparametersOptions desc("");
    desc.add(cmd_opts);

    if (parameterSet == "su3heatbath") {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12", "readMultipleConfs", "readFromConfNumber",
                                                       "readUntilConfNumber", "readConfsEvery"}))
            .add(ParametersIo::options.deleteSome({"fermObsInSingleFile", "fermObsCorrelatorsPrefix",
                                                   "fermObsCorrelatorsPostfix", "fermObsPbpPrefix", "fermObsPbpPostfix",
                                                   "hmcObsToSingleFile", "hmcObsPrefix", "hmcObsPostfix",
                                                   "rhmcObsToSingleFile", "rhmcObsPrefix", "rhmcObsPostfix"}))
            .add(ParametersGauge::options.keepOnlySome({"beta"}))
            .add(ParametersHeatbath::options)
            .add(ParametersObs::options.keepOnlySome({"measureTransportCoefficientKappa", "measureRectangles"}));
    } else if (parameterSet == "gaugeobservables") {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12"}))
            .add(ParametersIo::options.keepOnlySome(
                {"nDigitsInConfCheckpoint", "confPrefix", "confPostfix", "prngPrefix", "prngPostfix",
                 "rectanglesFilename", "transportCoefficientKappaFilename", "profilingDataPrefix",
                 "profilingDataPostfix", "gaugeObsInSingleFile", "gaugeObsPrefix", "gaugeObsPostfix"}))
            .add(ParametersGauge::options.keepOnlySome({"beta"}))
            .add(ParametersObs::options.keepOnlySome({"measureTransportCoefficientKappa", "measureRectangles"}));
    } else if (parameterSet == "inverter") {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12"}))
            .add(ParametersIo::options.deleteSome({"onlineMeasureEvery", "createCheckpointEvery",
                                                   "overwriteTemporaryCheckpointEvery", "hmcObsToSingleFile",
                                                   "hmcObsPrefix", "hmcObsPostfix", "rhmcObsToSingleFile",
                                                   "rhmcObsPrefix", "rhmcObsPostfix"}))
            .add(ParametersGauge::options.deleteSome({"gaugeAction"}))
            .add(ParametersFermion::options.deleteSome({"csw", "kappaMP", "muMP", "cswMP", "fermionActionMP"}))
            .add(ParametersSolver::options.deleteSome({"solverMP", "cgMaxIterationsMP", "restartEveryMP"}))
            .add(ParametersSources::options)
            .add(ParametersObs::options);
    } else if (parameterSet == "hmc") {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12", "readMultipleConfs", "readFromConfNumber",
                                                       "readUntilConfNumber", "readConfsEvery"}))
            .add(ParametersIo::options.deleteSome({"rhmcObsToSingleFile", "rhmcObsPrefix", "rhmcObsPostfix"}))
            .add(ParametersGauge::options.keepOnlySome({"beta", "gaugeAction"}))
            .add(ParametersFermion::options.deleteSome({"mass", "csw", "cswMP"}))
            .add(ParametersSolver::options)
            .add(ParametersSources::options)
            .add(ParametersObs::options)
            .add(ParametersHeatbath::options.keepOnlySome({"nThermalizationSteps"}))
            .add(ParametersIntegrator::options)
            .add(ParametersHmc::options);
    } else if (parameterSet == "rhmc") {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12", "readMultipleConfs", "readFromConfNumber",
                                                       "readUntilConfNumber", "readConfsEvery"}))
            .add(ParametersIo::options.deleteSome({"hmcObsToSingleFile", "hmcObsPrefix", "hmcObsPostfix"}))
            .add(ParametersGauge::options.keepOnlySome({"beta", "gaugeAction"}))
            .add(ParametersFermion::options.deleteSome(
                {"fermionActionMP", "kappa", "mu", "csw", "kappaMP", "muMP", "cswMP"}))
            .add(ParametersSolver::options.deleteSome(
                {"solverMP", "cgMaxIterationsMP", "restartEvery", "restartEveryMP"}))
            .add(ParametersSources::options.deleteSome({"placeSourcesOnHost"}))
            .add(ParametersObs::options)
            .add(ParametersIntegrator::options)
            .add(ParametersHeatbath::options.keepOnlySome({"nThermalizationSteps"}))
            .add(ParametersRationalApproximation::options)
            .add(ParametersRhmc::options);
    } else  // default: add all options
    {
        desc.add(ParametersConfig::options.deleteSome({"useReconstruct12"}))
            .add(ParametersIo::options)
            .add(ParametersGauge::options)
            .add(ParametersHeatbath::options)
            .add(ParametersFermion::options)
            .add(ParametersSolver::options)
            .add(ParametersSources::options)
            .add(ParametersObs::options)
            .add(ParametersHmc::options)
            .add(ParametersIntegrator::options)
            .add(ParametersRationalApproximation::options)
            .add(ParametersRhmc::options)
            .add(ParametersTest::options);
    }

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_opts).run(), vm);

    if (vm.count("help")) {  // See https://stackoverflow.com/a/5517755 as to why this is done before po::notifiy(vm)
        desc.printOptionsInCustomizedWay(std::cout);
        throw Inputparameters::help_required();
    }

    if (vm.count("inputFile")) {
        std::string config_file = vm["inputFile"].as<std::string>();
        ConfigFileNormalizer normalizer;
        add_option_aliases(&normalizer);
        // add stuff from input file
        std::string normalized_file;
        try {
            normalized_file = normalizer(config_file);
        } catch (std::invalid_argument) {
            std::cout << "Invalid config file " << config_file << std::endl;
            throw Inputparameters::parse_aborted();
        }
        std::istringstream normalized_file_stream(normalized_file);
        // todo: do not pass the "help" option...
        po::store(po::parse_config_file(normalized_file_stream, desc, false), vm);
    }

    po::notify(vm);  // checks whether all required arguments are set
    ChecksStringOptionsAndMapToEnum();
}

static void add_option_aliases(meta::ConfigFileNormalizer* const normalizer)
{
    normalizer->add_alias("NS", "nSpace");
    normalizer->add_alias("NT", "nTime");
    normalizer->add_alias("use_evenodd", "useEO");
    normalizer->add_alias("thermsteps", "nThermalizationSteps");
    normalizer->add_alias("thermalizationsteps", "nThermalizationSteps");
    normalizer->add_alias("puregauge", "useGaugeOnly");
    normalizer->add_alias("PGT", "useGaugeOnly");
    normalizer->add_alias("test_ref_value", "testRefVal");
    normalizer->add_alias("test_ref_value2", "testRefVal2");
    normalizer->add_alias("ThetaT", "thetaFermionTemporal");
    normalizer->add_alias("ThetaS", "thetaFermionSpatial");
}
