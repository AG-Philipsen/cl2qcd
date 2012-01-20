#!/usr/bin/env python
# coding=utf8
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

import csv
import optparse

if __name__ == '__main__':

	parser = optparse.OptionParser(usage='%prog [options] FILES...')
	parser.add_option('-p', '--plot', action='store_true', default=False, help='Plot the performance data')

	(args, files) = parser.parse_args()

	if len(files) < 1:
		files = ['inverter_benchmarks_profiling_data']


	lines = []
	started = False

	for filename in files:
		for line in open(filename, 'r'):
			if line.find('#device 0') != -1:
				started = True
			elif line.strip()[0] == '#':
				started = False
			elif started:
				lines.append(line)

	parsed = []
	for line in csv.reader(lines, delimiter='\t', skipinitialspace=True):
		parsed.append(line)

	sortedLines = sorted(parsed, key=lambda line: int(line[1]))
	totalTime = sum(map(lambda line: int(line[1]), sortedLines))

	if args.plot:
		import matplotlib.pyplot as plt
		import numpy as np

		tics = map(lambda line: line[0], sortedLines)
		ticpos = np.arange(len(tics)) + .5
		bars = map(lambda line: int(line[1]), sortedLines)

		rects = plt.barh(ticpos, bars)
		plt.yticks(ticpos, tics)
		plt.xlabel('Total Time (mus)')

		plt.show()

	else:
		print '{0:40}   {1:>10}   {2:>10}'.format('Kernel', 'Total Time [mus]', 'Rel. Time [mus]')
		for line in sortedLines:
			print '{0:40}   {1:>10}   {2:>10.1f}'.format(line[0], int(line[1]), int(line[1]) * 100. / totalTime)
