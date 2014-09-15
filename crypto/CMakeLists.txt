# Copyright 2013 Matthias Bach
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

find_package( Nettle REQUIRED )
include_directories(SYSTEM ${Nettle_INCLUDE_DIR})

# Build one crypto lib
add_library(crypto
	md5.cpp
)
target_link_libraries(crypto
	${Nettle_LIBRARIES}
)

add_executable(cryptoTest test.cpp)
target_link_libraries(cryptoTest ${Boost_LIBRARIES} crypto)
add_test(Crypto cryptoTest)
