#!/usr/bin/env python3
#
# Copyright (c) 2013 Christopher Pinke
# Copyright (c) 2014,2021 Alessandro Sciarra
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

import argparse
import os
import errno
import subprocess
import re

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Benchmark a kernel performance')
	parser.add_argument('-k', '--kernel', required=True, choices=['dks', 'dslash', 'su3heatbath'], help='Benchmark kernel, required')
	parser.add_argument('-p', '--path', default='./', help='Path to benchmark executable (default: ./')
	parser.add_argument('-n', '--numDevices', default=1, type=int, help='Number of devices to use (default: 1)')
	parser.add_argument('--ntMax', default=128, type=int, help='Maximal lattice extent in temporal direction (default: 128)')
	parser.add_argument('--ntMin', default=4, type=int, help='Minimal lattice extent in temporal direction (default: 4)')
	parser.add_argument('--ntIncr', default=4, type=int, help='Increment of lattice extent in temporal direction (default: 4)')
	parser.add_argument('--nsMax', default=64, type=int, help='Maximal lattice extent in spatial direction (default: 64)')
	parser.add_argument('--nsMin', default=16, type=int, help='Minimal lattice extent in spatial direction (default: 16)')
	parser.add_argument('--nsIncr', default=8, type=int, help='Increment of lattice extent in spatial direction (default: 8)')
	parser.add_argument('--calls', default=2000, type=int, help='Number of kernel calls (default: 2000)')
	parser.add_argument('--outFilePrefix', default="benchmark_", help='Prefix for output file (default: benchmark_)')
	parser.add_argument('--doNotUseGPU', default=False, action='store_true', help='Use CPU instead of GPU')
	parser.add_argument('-q', '--quite', default=False, action='store_true', help='Print benchmark output to standard output')
	args = parser.parse_args()

	# further setup based on command line options
	if args.kernel == 'dks':
		executable = args.path + 'dks_benchmark'
		labelOutput = 'D_KS'
	elif args.kernel == 'dslash':
		executable = args.path + 'dslash_benchmark'
		labelOutput = 'Dslash'
	elif args.kernel == 'su3heatbath':
		executable = args.path + 'su3heatbath_benchmarks'
		labelOutput = 'SU3_Heatbath'
	if args.doNotUseGPU == True:
		useGPU = 'false'
		useCPU = 'true'
	else:
		useGPU = 'true'
		useCPU = 'false'

	# open file for results
	f = open(f'{args.outFilePrefix}{args.kernel}_{args.numDevices}_devices.dat', "w")
	f.write("#Ns\tNt\tGFLOPS\t\tGB/s\t\ttime{msec}\n")

	for Ns in range(args.nsMin, args.nsMax+1, args.nsIncr):
		for Nt in range(args.ntMin, args.ntMax+1, args.ntIncr):
			# do not simulate all sizes...
			if Nt > 32:
				if not Nt % 32 == 0:
					continue
			print(f'\n# Testing {labelOutput} kernel on {Ns}^3*{Nt} lattice with {args.numDevices} device(s)..."')
			if not args.quite:
				print("# CL2QCD output:\n")

			# call program
			if not os.path.isfile(executable):
				raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), executable)

			p = subprocess.Popen([executable, f'--nSpace={Ns}', f'--nTime={Nt}', f'--useGPU={useGPU}', f'--useCPU={useCPU}', f'--nDevices={args.numDevices}', f'--nBenchmarkIterations={args.calls}'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			out, err = p.communicate()
			out = out.decode('utf8')
			lines = out.split("\n")
			err = err.decode('utf8')
			if not err == "":
				print("# There were errors:\n")
				print(err)
				if not args.quite:
					for line in lines:
						print(line)
				continue

			# get result from output
			performance = float('nan')
			memory = float('nan')
			time = float('nan')
			for line in lines:
				if not args.quite:
					print(line)
				if performance != performance:
					match = re.search(fr'{labelOutput} performance:\s+\d+.\d+', line)
					if match:
						performance = float(match.group().split()[2] )
				if memory != memory:
					match2 = re.search(fr'{labelOutput} memory:\s+\d+.\d+', line)
					if match2:
						memory = float(match2.group().split()[2] )
				if time != time:
					if args.numDevices == 1:
						match3 = re.search(fr'{labelOutput} device time per call:\s+\d+.\d+', line)
					else:
						match3 = re.search(fr'{labelOutput} host time per call:\s+\d+.\d+', line)
					if match3:
						time = float(match3.group().split()[5] )

			# Report
			print(f'# Measured performance of {performance} GFLOPS')
			print(f'# Measured memory of {memory} GB/S')
			if args.numDevices == 1:
				print(f'# Measured device time per call of {time} ms')
			else:
				print(f'# Measured hosttime per call of {time} ms')
			f.write(f'{Ns}\t{Nt}\t{performance}\t\t{memory}\t\t{time}\n')

	print('')
	f.close()
