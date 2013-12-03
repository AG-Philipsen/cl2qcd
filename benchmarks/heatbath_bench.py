#!/usr/bin/env python
#
# Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
# Matthias Bach, Christian Schaefer, Stefano Lottini, Alessandro Sciarra
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

from subprocess import *
import argparse
import os
import sys
from tempfile import NamedTemporaryFile

# benchmark configuration

EXECUTABLE = "heatbath_benchmarks"

NS = (16, 24, 32, 48)
NT = (4, 8, 12, 16, 24, 32, 48)

INPUT_TEMPLATE = """
##################################################################
##	This is 
##		heatbach_bench_input_{0}_{1}
##	It should reproduce the data obtained in 
##		hep-lat/9602007
##	There, NT=16, NS=16, beta=6.0 is reported to give
##		P: 0.593678 (24)
##	after 5000-10000 steps on the thermalized system.
##	Thus, the error in this setup is expected to be higher!
##	The output is written to the file
##		"gaugeobservables_beta6"
##################################################################
NS={0}
NT={1}
use_gpu={2}
beta=6
startcondition=cold
sourcefile=conf.save
thermalization=0
heatbathsteps=100
overrelaxsteps=1
savefrequency=10
print_to_screen=1
"""

# implementation

def generate_lattice_sizes():
	"""Generates combinations of NS and NT"""
	return ((ns, nt) for ns in NS for nt in NT)

def create_config_file(ns, nt, cpu=False):
	"""Create a configuration file for the given lattice size. If cpu flag ist not given
	   The configuration will allow GPU usage."""
	f = NamedTemporaryFile(delete=False);
	f.write(INPUT_TEMPLATE.format(ns, nt, 0 if cpu else 1))
	f.close()
	return f


def run_heatbath(input_file, verbose=False):
	"""Run the binary performing the heatbath.
	   If verbose is not set output of the executable will be suppressed."""
	cmd = ['./' + EXECUTABLE] + [input_file.name]

	if(verbose):
		subject =  Popen(cmd)
	else:
		subject =  Popen(cmd, stdout=PIPE)
		for line in subject.stdout:
			pass

	subject.wait()

	if subject.returncode == 0:
		print "BENCHMARK RUNNER: Program completed successfully"
	else:
		print "BENCHMARK RUNNER: Program terminated with exit code %i" % (subject.returncode )


def main():
	parser = argparse.ArgumentParser(description='Heatbath benchmark runner')
	parser.add_argument('--cpu', action='store_true', default=False, help='Run on CPU instead of GPU')
	parser.add_argument('--verbose', action='store_true', default=False, help='Show output of heatbath executable')
	args = parser.parse_args()

	for ns, nt in generate_lattice_sizes():
		print 'BENCHMARK RUNNER: Benchmarking {0}^3 * {1} lattice'.format(ns, nt)

		input_file = create_config_file(ns, nt, args.cpu)
		try:
			run_heatbath(input_file, args.verbose)
		finally:
			os.remove(input_file.name)


if __name__ == '__main__':
	sys.exit(main())
