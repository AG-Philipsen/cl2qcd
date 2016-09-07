/** @file
 * sourcesInterface.hpp
 *
 * Copyright 2016 Alessandro Sciarra
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

#pragma once

namespace physics {

    class SourcesParametersInterface {
        public:
            virtual ~SourcesParametersInterface(){}
            virtual bool placeSourcesOnHost() const = 0;
            virtual common::sourcetypes getSourceType() const = 0;
            virtual unsigned getSourceT() const = 0;
            virtual unsigned getSourceX() const = 0;
            virtual unsigned getSourceY() const = 0;
            virtual unsigned getSourceZ() const = 0;
            virtual unsigned getNt() const = 0;
            virtual unsigned getNs() const = 0;
    };

}
