#!/usr/bin/env python
#
# Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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
import re
import sys
import os

tolerance = 0.5

def compareFloat(reference, actual, precision):
	return actual > reference * (1.0 - precision / 100) and actual < reference * (1.0 + precision / 100)


def main():
	# Expect first arg to be reference file, remaining args are passed to application

	reference = open(sys.argv[1])
	#print '%f' % float(reference.readline())
	#print '%f' % float(reference.readline())

	#rm existing rhmc_output file
	# NOTE: the "-f" will suppress an error if the file does not exist
	os.system('rm -f rhmc_output')

	# perform RHMC with given input file
	# open tmp file to save the output of the rhmc
	# NOTE: stdout=PIPE does not work here, apparently no output file is created then
	subject = Popen(['../../rhmc'] + sys.argv[2:], stdout = PIPE)
	

	for line in subject.stdout:
		# Echo line to allow checking what's going on
		print line,

	subject.wait()
	if subject.returncode == 0 or subject.returncode == -11: # yes, it's kind of dirty to ignore -11 like that
		#print "Program completed successfully"
		pass
	else:
		print "Program terminated with exit code %i" % ( subject.returncode )
		return subject.returncode

	#now the value of interest is in the 1st row of "rhmc_output"
	candidate = open('rhmc_output')

	#get reference value, this is given in the first line of the reference file
	refval = float(reference.readline())
	val = float(candidate.readline().split()[1])
	candidate.close()
	reference.close()

	if not compareFloat(val, refval, tolerance):
		print "Invalid Result! Expected: %.16f  Real: %.16f (at %.1e %% tolerance)" % (refval, val, tolerance)
		return 127

	# no check failed
	print "Test done! No errors detected!"
	return 0

if __name__ == '__main__':
	sys.exit(main())
