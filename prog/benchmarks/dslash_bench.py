#!/usr/bin/env python

from subprocess import *
import re
import sys
import os
import math
import shutil
import dslash_bench_defs

def main():

	# some needed vars that are defined in the file dslash_bench_defs.py for convenience
	input_glob = dslash_bench_defs.input_glob
	input_var1 = dslash_bench_defs.input_var1
	input_var2 = dslash_bench_defs.input_var2
	debug = dslash_bench_defs.debug
	stdout = dslash_bench_defs.stdout
	backup = dslash_bench_defs.backup
	executable = dslash_bench_defs.executable

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
				f = open('tmpaaaaaaaa', 'w')
				f.write(input_glob + input_var1[iteration1] + input_var2[iteration2])
				f.close()
				args = ['./' + executable] + ['tmpaaaaaaaa']

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
					shutil.copy('tmpaaaaaaaa', backupname)

			if subject.returncode == 0:
				print "\tProgram completed successfully"
			else:
				print "\tProgram terminated with exit code %i" % (subject.returncode )
				continue

	if(switch == 0):
		os.remove('tmpaaaaaaaa')


if __name__ == '__main__':
	sys.exit(main())
