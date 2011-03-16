# Script to locate the LIME library
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
# Copyright 2011 Matthias Bach <bach@compeng.uni-frankfurt.de>

FIND_PACKAGE( PackageHandleStandardArgs )

IF (APPLE)

	FIND_LIBRARY(LIME_LIBRARIES lime DOC "LIME lib for OSX"
		ENV DYLD_LIBRARY_PATH
		ENV LIBRARY_PATH )
ELSE (APPLE)

	IF (WIN32)

		FIND_LIBRARY(LIME_LIBRARIES lime
			ENV PATH
		)

	ELSE (WIN32)

		# Unix style platforms
		FIND_LIBRARY(LIME_LIBRARIES lime
			ENV LD_LIBRARY_PATH
			ENV LIBRARY_PATH
		)

	ENDIF (WIN32)

ENDIF (APPLE)

FIND_PATH(LIME_INCLUDE_DIR lime.h DOC "Include for LIME")

# Also search relative to lib ( git build )
IF ( NOT LIME_INCLUDE_DIR )
	GET_FILENAME_COMPONENT(_LIME_LIB_DIR ${LIME_LIBRARIES} PATH)
	GET_FILENAME_COMPONENT(_LIME_INC_CAND ${_LIME_LIB_DIR}/../include ABSOLUTE)
	FIND_PATH(LIME_INCLUDE_DIR lime.h PATHS ${_LIME_INC_CAND} )
ENDIF ( NOT LIME_INCLUDE_DIR )

FIND_PACKAGE_HANDLE_STANDARD_ARGS( LIME DEFAULT_MSG LIME_LIBRARIES LIME_INCLUDE_DIR )

