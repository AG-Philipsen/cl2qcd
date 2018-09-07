# Build specifications for the contents of the hardware package
#
# Copyright (c) 2015 Christopher Pinke
# Copyright (c) 2016,2018 Alessandro Sciarra
#
# This file is part of CL2QCD.
#
# CL2QCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CL2QCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.

add_unit_test(NAME interfacesImplementations/hardwareParameters      LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/openClKernelParameters  LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/physicsParameters       LIBRARIES meta exceptions)
add_unit_test(NAME interfacesImplementations/observablesParameters   LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/latticesParameters      LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/fermionmatrixParameters LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/algorithmsParameters    LIBRARIES meta)
add_unit_test(NAME interfacesImplementations/interfacesHandler       LIBRARIES meta exceptions)
