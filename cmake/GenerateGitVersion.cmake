# Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
# Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.

find_package(Git)

if(GIT_FOUND)
	execute_process( COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:%H
	                 WORKING_DIRECTORY ${SOURCE_DIR}
	                 OUTPUT_VARIABLE GIT_COMMIT_ID
	                 OUTPUT_STRIP_TRAILING_WHITESPACE )
	execute_process( COMMAND ${GIT_EXECUTABLE} diff --shortstat
	                 WORKING_DIRECTORY ${SOURCE_DIR}
	                 OUTPUT_VARIABLE GIT_DIRTY
	                 OUTPUT_STRIP_TRAILING_WHITESPACE )
	if( GIT_DIRTY )
		set( GIT_COMMIT_ID "${GIT_COMMIT_ID}*" )
	endif( GIT_DIRTY )
	message("Git commit id: ${GIT_COMMIT_ID}")
else(GIT_FOUND)
	message(WARNING "Failed to find git, can't write git commit id to executables")
	set(GIT_COMMIT_ID "N/A")
endif(GIT_FOUND)
file(WRITE gitcommitid.h.new "#define GIT_COMMIT_ID \"${GIT_COMMIT_ID}\"\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
		gitcommitid.h.new gitcommitid.h)
