# Build configuration for OpenCL based HMC
#
# Copyright 2011 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
# Copyright 2011-2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
# Copyright 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

# the folowing are pure profiling benchmarks
add_executable(heatbath_benchmarks heatbath_benchmarks.cpp)
target_link_libraries(heatbath_benchmarks optimal)

add_executable(dslash_benchmark dslash_benchmark.cpp)
target_link_libraries(dslash_benchmark optimal)

add_executable(dslash_benchmark2 dslashBenchmarkMain.cpp)
target_link_libraries(dslash_benchmark2 optimal)

add_executable(dslash_multidev dslash_multidev.cpp)
target_link_libraries(dslash_multidev optimal)

add_executable(dks_multidev dks_multidev.cpp)
target_link_libraries(dks_multidev optimal)
