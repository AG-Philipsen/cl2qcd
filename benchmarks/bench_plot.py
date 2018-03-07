# coding=utf8
#
# Copyright (c) 2012,2013 Matthias Bach
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

import matplotlib.pyplot as plt
import numpy as np
import re
from collections import namedtuple

linestyles_color = ['r.', 'b.', 'r*', 'b*', 'g.', 'k.', 'r,', 'b,', 'g,', 'k,', 'g*', 'k*']
linestyles_gray = ['kd', 'kh', 'ko', 'k*', 'k^', 'kv', 'kx', 'k+']
#linestyles_gray = ['k+', 'kx', 'k1', 'k2', 'k3', 'k4', 'kd', 'k.']

FileData = namedtuple('FileData', ['label', 'runs', 'xpos'])

def main(datafiles, filelabels, kernelpattern, output=None, metric='both', title=False, maxSize=None, legend_pos = None, offset_lines = False, grayscale=False):

	linestyles = linestyles_gray if grayscale else linestyles_color;
	markersize = 10 if grayscale else 15

	filedatas = []

	for i in xrange(len(datafiles)):
		runs = []

		for line in open(datafiles[i], 'r'):
			match = re.match(r'## NSPACE:\s+(\d+)', line)
			if match:
				nspace = int(match.group(1))
			match = re.match(r'## NTIME:\s+(\d+)', line)
			if match:
				ntime = int(match.group(1))
			floatpattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
			match = re.match(r'\s*(' + kernelpattern + ')\s+\d+\s+\d+\s+\d+\s+\d+\s+(?P<bandwidth>' + floatpattern + ')\s(?P<gflops>' + floatpattern + ')', line)
			if match:
				if maxSize and nspace**3 * ntime > maxSize:
					continue

				bandwidth = float(match.group('bandwidth'))
				gflops = float(match.group('gflops'))
				runs.append((nspace, ntime, bandwidth, gflops))

		# dump data to cli
		runs = np.array(sorted(runs, key=lambda p: p[0]**3 * p[1]))
		print 'NSPACE   NTIME    GB/S   GFLOPS'
		for run in runs:
			print '{0[0]:>5}   {0[1]:>5}   {0[2]:>5}   {0[3]:>5}'.format(run)

		xpos = map(lambda (x,y): int(x)**3 * int(y), runs[:,0:2])
		filedatas.append(FileData(filelabels[i], runs, xpos))

	xtic_pos = []
	for data in filedatas:
		xtic_pos.extend(map(lambda p, (x, y): (p, x, y), data.xpos, data.runs[:,0:2]))
	xtic_pos.sort(key=lambda e: e[0])
	print xtic_pos
	# make sure tics don't overlapp
	max_xtics = 30
	min_delta = (xtic_pos[-1][0] - xtic_pos[0][0]) / max_xtics
	i = len(xtic_pos) - 1
	while i > 0:
		if xtic_pos[i][0] < xtic_pos[i-1][0] + min_delta:
			del xtic_pos[i-1]
		i -= 1
	print xtic_pos
	xtic_label = map(lambda (p, x, y): '{0}^3x{1}'.format(int(x), int(y)), xtic_pos)
	latex_label = [r'\({0}^3 \times {1}\)'.format(int(x), int(y)) for (p,x,y) in xtic_pos]
	xtic_pos = map(lambda (p, x, y): p, xtic_pos)
	print 'xtic={{{0}}},'.format(','.join((str(tic) for tic in xtic_pos)))
	print 'xticlabels={{{0}}},'.format(','.join(latex_label))

	fig = plt.figure(figsize=(10,4))
	fig.subplots_adjust(bottom=0.28)
	ax1 = fig.add_subplot(111)
	if metric == 'both':
		ax2 = ax1.twinx()
	lines = []
	labels = []
	linestyle = 0

	if offset_lines:
		offset_per_line = (min_delta / 2)
		offset = - int(len(filedatas) / 2) * offset_per_line / 2
	else:
		offset_per_line = 0
		offset = 0

	for data in filedatas:
		xpos = [pos + offset for pos in data.xpos]
		if metric == 'gflops':
			lines.append(ax1.plot(xpos, data.runs[:,3], linestyles[linestyle], markersize=markersize))
			labels.append(data.label)
			linestyle += 1
		elif metric == 'gbytes':
			lines.append(ax1.plot(xpos, data.runs[:,2], linestyles[linestyle], markersize=markersize))
			labels.append(data.label)
			linestyle += 1
		else:
			lines.append(ax1.plot(xpos, data.runs[:,2], linestyles[linestyle], markersize=markersize))
			labels.append(data.label + ' Bandwidth')
			linestyle += 1
			lines.append(ax2.plot(xpos, data.runs[:,3], linestyles[linestyle], markersize=markersize))
			labels.append(data.label + ' Gflops')
			linestyle += 1
		offset += offset_per_line

	if title:
		ax1.set_title(title)
	ax1.set_xticks(xtic_pos)
	ax1.set_xticklabels(xtic_label, rotation=90)
	ax1.set_xlabel('Lattice Size')
	if metric in ('both', 'gbytes'):
		ax1.set_ylabel('Bandwidth GB/s')
	if metric == ('both'):
		ax2.set_ylabel('Gflops')
	if metric == ('gflops'):
		ax1.set_ylabel('Gflops')
	ax1.set_ylim(bottom=1)
	if metric == 'both':
		max_lim = (0, max(ax1.get_ylim()[1], ax2.get_ylim()[1]))
		ax1.set_ylim(max_lim)
		ax2.set_ylim(max_lim)

	extra_legend_args = {}
	if legend_pos:
		extra_legend_args['loc'] = legend_pos
	legend = fig.legend(lines, labels, **extra_legend_args)

	if output:
		plt.savefig(output)
	else:
		plt.show()
