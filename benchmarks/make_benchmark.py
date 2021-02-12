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
	parser.add_argument('-k', '--kernel', required=True, choices=['dks', 'dslash', 'su3heatbath'], help='Kernel to benchmark, required.')
	parser.add_argument('-p', '--path', default='./', help='Path to benchmark executable (default: "./").')
	parser.add_argument('-s', '--nSpace', nargs='+', type=int, help='Spacial extend of the lattice to use. It accepts multiple values.')
	parser.add_argument('-t', '--nTime', nargs='+', type=int, help='Temporal extend of the lattice to use. It accepts multiple values.')
	parser.add_argument('--calls', default=2500, type=int, help='Number of kernel calls (default: 2500).')
	parser.add_argument('--outFilePrefix', default="benchmark_", help='Prefix for output file (default: "benchmark_").')
	parser.add_argument('--doNotUseGPU', default=False, action='store_true', help='Use CPU instead of GPU.')
	parser.add_argument('-q', '--quite', default=False, action='store_true', help='Print benchmark output to standard output.')
	parser.add_argument('--keepAuxFiles', default=False, action='store_true', help='Do not remove CL2QCD auxiliary files.')

	group = parser.add_mutually_exclusive_group()
	group.add_argument('-n', '--numDevices', default=1, type=int, help='Number of devices to use.')
	group.add_argument('-d', '--devices', nargs='+', type=int, help='Which devices to use.')

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

	space_dims = args.nSpace if args.nSpace else [16, 24, 32, 48]
	time_dims = args.nTime if args.nTime else [4, 8, 12, 16, 24, 32, 64, 96, 128]

	# open file for results
	f = open(f'{args.outFilePrefix}{args.kernel}_{args.numDevices}_devices.dat', "w")
	f.write("#{0:9}{1:15}{2:15}{3:15}{4:10}\n".format('Ns', 'Nt', 'GFLOPS', 'GB/s', 'ms/call'))

	cl2qcd_command_line_common_options = [executable, f'--useGPU={useGPU}', f'--useCPU={useCPU}', f'--nBenchmarkIterations={args.calls}']
	if args.numDevices != None:
		cl2qcd_command_line_common_options += [f'--nDevices={args.numDevices}']
	if args.devices != None:
		for id in args.devices:
			cl2qcd_command_line_common_options += [f'--deviceId={id}']

	benchmark_number = 0
	total_number_of_benchmarks = len(space_dims) * len(time_dims)
	for Ns in space_dims:
		for Nt in time_dims:
			benchmark_number += 1
			print(f'\n# Benchmark {benchmark_number} of {total_number_of_benchmarks}')
			print(f'# Testing {labelOutput} kernel on {Nt}x{Ns}^3 lattice on {args.numDevices} device(s) [{args.calls} calls]')
			if not args.quite:
				print("# CL2QCD output:\n")

			# call program
			if not os.path.isfile(executable):
				raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), executable)

			cl2qcd_command_line_options = cl2qcd_command_line_common_options + [f'--nSpace={Ns}', f'--nTime={Nt}']
			p = subprocess.Popen(cl2qcd_command_line_options, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
			f.write(f'{Ns:<10d}{Nt:<15}{performance:<15.3f}{memory:<15.3f}{time:<10.3f}\n')

	print('')
	f.close()

	# Remove auxiliary files of CL2QCD
	if not args.keepAuxFiles:
		try:
			os.remove(f'{args.kernel}_benchmark.log')
			os.remove("general_time_output")
			os.remove(f'{args.kernel}_benchmark_profiling_data')
		except OSError as e:  ## if failed, report it back to the user ##
			print (f'Error: {e.filename} - {e.strerror}.')
