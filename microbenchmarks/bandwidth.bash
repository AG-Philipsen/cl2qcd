#!/bin/bash

# Copyright (c) 2011 Matthias Bach
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

MIN_NATIVE_GROUP_SIZE=64
MAX_NATIVE_GROUP_SIZE=256

NUM_COMPUTE_UNITS=24

if (( $# < 1 )); then
	echo "ERROR: Missing required argument"
	echo "Usage bandwidth.sh <result_file_prefix>"
	exit -1
fi

RESULT_PREFIX=$1

./bandwidth --threads=$MIN_NATIVE_GROUP_SIZE | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.minthreads.float.chunk.dat"
./bandwidth --threads=$MAX_NATIVE_GROUP_SIZE | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.maxthreads.float.chunk.dat"

./bandwidth --threads=$MIN_NATIVE_GROUP_SIZE --su3 | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.minthreads.su3.chunk.dat"
./bandwidth --threads=$MAX_NATIVE_GROUP_SIZE --su3 | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.maxthreads.su3.chunk.dat"

./bandwidth --threads=$MIN_NATIVE_GROUP_SIZE --single | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.minthreads.float.single.dat"
./bandwidth --threads=$MAX_NATIVE_GROUP_SIZE --single | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.maxthreads.float.single.dat"

./bandwidth --threads=$MIN_NATIVE_GROUP_SIZE --su3 --single | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.minthreads.su3.single.dat"
./bandwidth --threads=$MAX_NATIVE_GROUP_SIZE --su3 --single | grep -v '^\[' > "$RESULT_PREFIX.sweepgroups.maxthreads.su3.single.dat"


./bandwidth --groups=$NUM_COMPUTE_UNITS --threads=$MAX_NATIVE_GROUP_SIZE | grep -v '^\[' > "$RESULT_PREFIX.sweepthreads.maxthreads.float.chunk.dat"
./bandwidth --groups=$NUM_COMPUTE_UNITS --threads=$MAX_NATIVE_GROUP_SIZE --su3 | grep -v '^\[' > "$RESULT_PREFIX.sweepthreads.maxthreads.su3.chunk.dat"
./bandwidth --groups=$NUM_COMPUTE_UNITS --threads=$MAX_NATIVE_GROUP_SIZE --single | grep -v '^\[' > "$RESULT_PREFIX.sweepthreads.maxthreads.float.single.dat"
./bandwidth --groups=$NUM_COMPUTE_UNITS --threads=$MAX_NATIVE_GROUP_SIZE --su3 --single | grep -v '^\[' > "$RESULT_PREFIX.sweepthreads.maxthreads.su3.single.dat"
