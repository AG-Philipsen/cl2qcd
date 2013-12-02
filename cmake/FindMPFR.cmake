# Script to locate the MPFR library
#
# This script is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This script is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this script.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>

FIND_PACKAGE(PackageHandleStandardArgs)

IF(APPLE)
	FIND_LIBRARY(MPFR_LIBRARIES mpfr DOC "MPFR lib for OSX"
		ENV DYLD_LIBRARY_PATH
		ENV LIBRARY_PATH )
ELSE()
	IF(WIN32)
		FIND_LIBRARY(MPFR_LIBRARIES mpfr
			ENV PATH
		)
	ELSE()
		# Unix style platforms
		FIND_LIBRARY(MPFR_LIBRARIES mpfr
			ENV LD_LIBRARY_PATH
			ENV LIBRARY_PATH
		)
	ENDIF()
ENDIF()

FIND_PATH(MPFR_INCLUDE_DIR mpfr.h DOC "Include for MPFR")

# Also search relative to lib ( git build )
IF ( NOT MPFR_INCLUDE_DIR )
	GET_FILENAME_COMPONENT(_MPFR_LIB_DIR ${MPFR_LIBRARIES} PATH)
	GET_FILENAME_COMPONENT(_MPFR_INC_CAND ${_MPFR_LIB_DIR}/../include ABSOLUTE)
	FIND_PATH(MPFR_INCLUDE_DIR mpfr.h PATHS ${_MPFR_INC_CAND} )
ENDIF ( NOT MPFR_INCLUDE_DIR )

FIND_PACKAGE_HANDLE_STANDARD_ARGS( MPFR DEFAULT_MSG MPFR_LIBRARIES MPFR_INCLUDE_DIR )
