#!/usr/bin/env python

from subprocess import *
import re
import sys
import os

tolerance = 0.1

def compareFloat(reference, actual, precision):
	return actual > reference * (1.0 - precision / 100) and actual < reference * (1.0 + precision / 100)


def main():
	# Expect first arg to be reference file, remaining args are passed to application

	reference = open(sys.argv[1])
	#print '%f' % float(reference.readline())
	#print '%f' % float(reference.readline())

	#rm existing hmc_output file
	# NOTE: the "-f" will suppress an error if the file does not exist
	os.system('rm -f hmc_output')

	# perform HMC with given input file
	# open tmp file to save the output of the hmc
	# NOTE: stdout=PIPE does not work here, apparently no output file is created then
	subject = Popen(['../hmc'] + sys.argv[2:], stdout = PIPE)

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

	#now the value of interest is in the 2nd row of "hmc_output"
	os.system('awk \'{print $2}\' hmc_output > hmc_test_tmp')
	candidate = open('./hmc_test_tmp', 'r')

	#get reference value, this is given in the first line of the reference file
	refval = float(reference.readline())
	val = float(candidate.readline())
	candidate.close()
	reference.close()

	if not compareFloat(val, refval, tolerance):
		print "Invalid Result! Expected: %f  Real: %f (at %.1f %% tolerance)" % (refval, val, tolerance)
		return 127

	# rm tmp files
        os.system('rm -f hmc_hmc_test_tmp')


	# no check failed
	return 0

if __name__ == '__main__':
	sys.exit(main())
