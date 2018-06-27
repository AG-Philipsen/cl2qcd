/** @file
 * Input file handling implementation
 *
 * Copyright (c) 2018 Francesca Cuteri
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

#include "parametersBasic.hpp"

#include "../executables/exceptions.hpp"

#include <sys/ioctl.h>

static unsigned short int getTerminalWidth()
{
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
}

unsigned short int getHelperWidth()
{
    unsigned short int threshold = getTerminalWidth() / 10 * 9;
    return (threshold > 210) ? 210 : threshold;
}

meta::InputparametersOptions::InputparametersOptions(std::string optionsDescription)
    : meta::InputparametersOptions(optionsDescription, getHelperWidth(), getHelperWidth() / 3)
{
}

meta::InputparametersOptions::InputparametersOptions(std::string optionsDescriptionIn, unsigned int lineLengthIn,
                                                     unsigned int minimumDescriptionLengthIn)
    : po::options_description(optionsDescriptionIn, lineLengthIn, minimumDescriptionLengthIn), lineLength(lineLengthIn)
{
}

void meta::InputparametersOptions::printOptionsInCustomizedWay(std::ostream& stream) const
{
    unsigned short int maxOptionWidth = get_option_column_width();
    unsigned short int optionWidth    = (maxOptionWidth < lineLength / 4) ? lineLength / 4 : maxOptionWidth + 5;
    print(stream, optionWidth);
}

meta::InputparametersOptions&
meta::InputparametersOptions::keepOnlySomeOptions(std::initializer_list<std::string> whichOptions)
{
    if (whichOptions.size() == 0)
        throw Print_Error_Message("No options to be kept passed!", __FILE__, __LINE__);

    std::vector<boost::shared_ptr<po::option_description>>&
        listOfOptions = const_cast<std::vector<boost::shared_ptr<po::option_description>>&>(options());
    std::vector<boost::shared_ptr<po::option_description>> newList{};
    for (auto requiredOption : whichOptions) {
        auto it = std::find_if(listOfOptions.begin(), listOfOptions.end(),
                               [&requiredOption](boost::shared_ptr<po::option_description> ptr) {
                                   return ptr->canonical_display_name() == requiredOption;
                               });
        if (it == listOfOptions.end())
            throw Print_Error_Message("Required option \"" + requiredOption + "\" not found!", __FILE__, __LINE__);
        else
            newList.push_back(*it);
    }
    listOfOptions = std::move(newList);
    return *this;
}
