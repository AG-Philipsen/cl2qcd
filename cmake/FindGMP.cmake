# Script to locate the GMP library
#
# Copyright (c) 2013 Matthias Bach
# Copyright (c) 2018 Alessandro Sciarra
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

FIND_PACKAGE(PackageHandleStandardArgs)

IF(APPLE)
	FIND_LIBRARY(GMP_LIBRARIES gmp DOC "GMP lib for OSX"
		ENV DYLD_LIBRARY_PATH
		ENV LIBRARY_PATH )
ELSE()
	IF(WIN32)
		FIND_LIBRARY(GMP_LIBRARIES gmp
			ENV PATH
		)
	ELSE()
		# Unix style platforms
		FIND_LIBRARY(GMP_LIBRARIES gmp
			ENV LD_LIBRARY_PATH
			ENV LIBRARY_PATH
		)
	ENDIF()
ENDIF()

FIND_PATH(GMP_INCLUDE_DIR gmp.h DOC "Include for GMP")

# Also search relative to lib ( git build )
IF ( NOT GMP_INCLUDE_DIR )
	GET_FILENAME_COMPONENT(_GMP_LIB_DIR ${GMP_LIBRARIES} PATH)
	GET_FILENAME_COMPONENT(_GMP_INC_CAND ${_GMP_LIB_DIR}/../include ABSOLUTE)
	FIND_PATH(GMP_INCLUDE_DIR gmp.h PATHS ${_GMP_INC_CAND} )
ENDIF ( NOT GMP_INCLUDE_DIR )

FIND_PACKAGE_HANDLE_STANDARD_ARGS( GMP DEFAULT_MSG GMP_LIBRARIES GMP_INCLUDE_DIR )
