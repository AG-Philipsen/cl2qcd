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

#ifndef _META_PARAMETERS_BASIC_HPP_
#define _META_PARAMETERS_BASIC_HPP_

#include "../common_header_files/types.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <cmath>
#include <iomanip>
#include <string>
#include <vector>
namespace po = boost::program_options;

namespace meta {

    class InputparametersOptions : public po::options_description {
      public:
        InputparametersOptions(std::string optionsDescription);
        void printOptionsInCustomizedWay(std::ostream& stream) const;
        /**
         * These methods are meant to be used in meta::Inputparameters in order to possibly extract only some options
         * from the parent ones, because in some executables not all the options of a parent are meaningful. To
         * implement these methods, one has to deal with the boost implementation of po::options_description which
         * cannot be changed. In particular, such a class contain some constant members and this automatically forbids
         * any copy/assignment of the class. We therefore use the po::options_description::options method to get a
         * reference to the options and we delete there those we do not want to keep. Unfortunately, this is not
         * straightforward, since the po::options_description::options method returns a constant reference and, hence,
         * we must cast away const.
         *
         * ATTENTION: Reading the po::options_description implementation and, in particular, the two methods
         *            po::options_description::add, it should be clear that an object of type po::options_description
         *            is more complicated than a simple collection of options, since it contains a private member of
         *            type std::vector<boost::shared_pointer<po::options_description>> and this is to allow to collect
         *            group of options (like we do in meta::Inputparameters). In general, the array containing all the
         *            single options is a union of all groups plus the options which are not in any group but which
         *            belong to the class itself. Since the array of options is indeed an object of type
         *            std::vector<boost::shared_pointer<po::option_description>> (note the singular in
         *            option_description), then when the group vector is nonempty, there will be pairs of shared
         *            pointers pointing to the same po::option_description (one in the group and one in the options).
         *            This breaks down the implementation of this method, which deletes only options from the options
         *            collection.
         *
         * CONCLUSION: Use it to delete options from the class collection ONLY, and not to delete options from collected
         *             groups! In meta::Inputparameters, for instance, this method is used on the parent's member of
         *             type InputparametersOptions, only.
         */
        InputparametersOptions& keepOnlySome(std::initializer_list<std::string> whichOptions);
        InputparametersOptions& deleteSome(std::initializer_list<std::string> whichOptions);

      private:
        InputparametersOptions(std::string optionsDescriptionIn, unsigned int lineLengthIn,
                               unsigned int minimumDescriptionLengthIn);
        unsigned int lineLength;  // Store even if privately present in parent
    };

    template<typename T>
    std::string getDefaultForHelper(const T& number)
    {
        // See https://stackoverflow.com/a/1736040 to understand for what this function is.
        std::ostringstream stream;
        if (std::abs(number) <= 1.e-5 || std::abs(number) >= 1.e+5) {
            stream << number;
            return stream.str();
        } else {
            stream << std::fixed << std::setprecision(16) << number;  // stream has always 16 decimal digits!
            std::string toBePrinted = stream.str();
            boost::trim_right_if(toBePrinted, boost::is_any_of("0"));
            boost::trim_right_if(toBePrinted, boost::is_any_of("."));
            return toBePrinted;
        }
    }

}  // namespace meta

#endif
