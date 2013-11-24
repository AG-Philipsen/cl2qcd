#!/usr/bin/env python
#
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
import re
import sys

tolerance = 0.1

def compareFloat(reference, actual, precision):
	return actual > reference * (1.0 - precision / 100) and actual < reference * (1.0 + precision / 100)


def main():
	# Expect first arg to be reference file, remaining args are passed to application

	reference = open(sys.argv[1])
	#print '%f' % float(reference.readline())
	#print '%f' % float(reference.readline())

	subject = Popen(['../inverter'] + sys.argv[2:], stdout=PIPE)

	for line in subject.stdout:
		# Echo line to allow checking what's going on
		print line,

		# Now we can check what the line looks like and react on it
		# The regexp searches for "[<timestamp>] INFO: <digits><whitespace><float>" and returns only the <float>-part.
		# For float the scientific notation matcher from http://docs.python.org/library/re.html#re.search is used.
		match = re.match(r"^\[\d\d:\d\d:\d\d\] INFO: \d+\s+([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$",line)
		if match != None:
			matched = float(match.group(1))
			#print 'Matched: %f' % matched
			try:
				refval = float(reference.readline())
			except ValueError:
				print "Too many values in output."
				return 126
			if not compareFloat(matched, refval, tolerance):
				print "Invalid Result! Expected: %f  Real: %f (at %.1f %% tolerance)" % (refval, matched, tolerance)
				return 127

	#print '--'

	subject.wait()
	if subject.returncode == 0:
		#print "Program completed successfully"
		pass
	else:
		print "Program terminated with exit code %i" % ( subject.returncode )
		return subject.returncode

	# check whether reference has been used up
	if reference.readline() != '':
		print "Not enough values have been reproduced."
		return 124

	# no check failed
	return 0

if __name__ == '__main__':
	sys.exit(main())
