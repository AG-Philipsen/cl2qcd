# Script to locate the Nettle library
#
# Copyright (c) 2013 Matthias Bach
# Copyright (c) 2017,2018 Alessandro Sciarra
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

FIND_PACKAGE( PackageHandleStandardArgs )

IF (APPLE)

	FIND_LIBRARY(Nettle_LIBRARIES nettle DOC "Nettle lib for OSX"
		ENV DYLD_LIBRARY_PATH
		ENV LIBRARY_PATH )
ELSE (APPLE)

	IF (WIN32)

		FIND_LIBRARY(Nettle_LIBRARIES nettle
			ENV PATH
		)

	ELSE (WIN32)

		# Unix style platforms
		FIND_LIBRARY(Nettle_LIBRARIES nettle
			ENV LD_LIBRARY_PATH
			ENV LIBRARY_PATH
		)

	ENDIF (WIN32)

ENDIF (APPLE)

FIND_PATH(Nettle_INCLUDE_DIR nettle/nettle-meta.h DOC "Include for Nettle")

# Also search relative to lib ( git build )
IF ( NOT Nettle_INCLUDE_DIR )
        GET_FILENAME_COMPONENT(_Nettle_LIB_DIR ${Nettle_LIBRARIES} PATH)
        GET_FILENAME_COMPONENT(_Nettle_INC_CAND ${_Nettle_LIB_DIR}/../include ABSOLUTE)
        FIND_PATH(Nettle_INCLUDE_DIR nettle/nettle-meta.h PATHS ${_Nettle_INC_CAND} )
ENDIF ( NOT Nettle_INCLUDE_DIR )

FIND_PACKAGE_HANDLE_STANDARD_ARGS( Nettle DEFAULT_MSG Nettle_LIBRARIES Nettle_INCLUDE_DIR )
