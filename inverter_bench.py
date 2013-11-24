#!/usr/bin/env python

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

from subprocess import *
import sys
import os
import shutil
import inverter_bench_defs
from tempfile import NamedTemporaryFile

def main():

	# some needed vars that are defined in the file inverter_bench_defs.py for convenience
	input_glob = inverter_bench_defs.input_glob
	input_var1 = inverter_bench_defs.input_var1
	input_var2 = inverter_bench_defs.input_var2
	debug = inverter_bench_defs.debug
	stdout = inverter_bench_defs.stdout
	backup = inverter_bench_defs.backup
	executable = inverter_bench_defs.executable

	switch = 0
	# check if input-file is given at argv[1]
	if(len(sys.argv) != 2):
		print "no input- or reference-file given. Perform all existing benchmarks!"
	else:
		input_name = sys.argv[1]
		print "\tbenchmarking \"" + executable + "\" using inputfile \"" + sys.argv[1]
		switch = 1

	# performs as many tests as specified in the defs-file
	if(switch == 0):
		size1 = len(input_var1)
		size2 = len(input_var2)
	else :
		size1 = 1
		size2 = 1

	for iteration1 in range(size1):
		for iteration2 in range(size2):
			iteration = iteration1*size2 + iteration2
			size = size1*size2
			print '\tbenchmark %i of %i'% (iteration+1,  size)

			# select input-file if wished, otherwise create one
			if(switch == 1):
				inputfile = sys.argv[1:]
				args = ['./' + executable] + sys.argv[1:]
			else:
				# open file "tmpaaaaaaaa", deleting whatever content was in it before
				f = NamedTemporaryFile(delete=False);
				f.write(input_glob + input_var1[iteration1] + input_var2[iteration2])
				f.close()
				args = ['./' + executable] + [f.name]

			# run the prog
			if(stdout):
				subject =  Popen(args)
			else:
				subject =  Popen(args, stdout=PIPE)

			subject.wait()
			if(switch == 0):
				if(backup):
					# save inputfile to a different file
					backupname = 'input_' + str(iteration)
					shutil.copy(f.name, backupname)

			if subject.returncode == 0:
				print "\tProgram completed successfully"
			else:
				print "\tProgram terminated with exit code %i" % (subject.returncode )
				continue

	if(switch == 0):
		os.remove(f.name)


if __name__ == '__main__':
	sys.exit(main())

