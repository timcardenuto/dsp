#! /usr/bin/env python

#------------------------------------------------------------------------------
# MIT License
#
# Copyright (C) 2017  Tim Cardenuto
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# For questions contact the originator at timcardenuto@gmail.com
#------------------------------------------------------------------------------

from __future__ import print_function
import numpy as np
import sys


# taps are vector of tap locations (bit/index locations)
# Ex. the polynomial x4 + x + 1 has taps [1 0]
# Ex. the polynomial x4 + x3 + 1 has taps [3 0]
# the order of the taps doesn't matter - they all get xor'd together (except the first one, the order)
def generic(order, taps, fill):
	N = pow(2,order) - 1;	# num of codes, also sequence repeat length in bits

	# chipping sequence comes out the lefts side of LFSR (zero bit)
	seq = np.zeros((N), dtype = np.int)
	newfill = np.zeros((order), dtype = np.int)

	# loop (aka clock) for one sequence length
	for i in range(N):
		seq[i] = fill[0]

		# shift all bits based on LFSR order
		for j in range(order-1):
			newfill[j] = fill[j+1]
		# add the 'new' bit at the end based on taps
		newfill[order-1] = 0
		for tap in taps:
			newfill[order-1] = newfill[order-1] ^ fill[tap]

		np.copyto(fill,newfill)

	return seq



def four(fill):
	if (len(fill) != 4):
		print("Error: this LFSR only accepts 4 bit fills")
		sys.exit(1)

	L = 4;          	# num of bits, or shift registers
	N = pow(2,L) - 1;	# num of codes, also sequence repeat length in bits

	# chipping sequence comes out the lefts side of LFSR (zero bit)
	seq = np.zeros((N), dtype = np.int)

	# loop (aka clock) for one sequence length
	for i in range(N):
		seq[i] = fill[0]
		# 1 + x3 + x4
		#newfill = [fill[1], fill[2], fill[3], (fill[3] ^ fill[0])]
		# 1 + x + x4
		newfill = [fill[1], fill[2], fill[3], (fill[1] ^ fill[0])]
		fill = newfill

	return seq


def generateFills(order):
	length = pow(2,order)
	fills = np.zeros((length,order), dtype = np.int)
	# loop through columns for length order
	# loop through alternating sequence in each column based on order

	for i in range(order):
		# this will treat array as backwards bytes, index 0 is LSB
		repeatsize = pow(2,i) # the bit position 'value' gives you the number of repeating values, for ex, the 4th bit is '8' so there are 8 zeros and then 8 ones
		count = 0
		state = 0
		for j in range(length):
			if (count < repeatsize and state == 0):
				fills[j][i] = 0
				count = count + 1
			elif (count < repeatsize and state == 1):
				fills[j][i] = 1
				count = count + 1
			else:
				count = 1
				if (state == 0):
					fills[j][i] = 1
					state = 1
				else:
					fills[j][i] = 0
					state = 0

	# delete the all zeros row, it's never used in LFSR's
	fills = np.delete(fills, (0), axis=0)
	# flip array to read like you'd expect, not reversed
	fills = np.fliplr(fills)
	return fills


if __name__ == "__main__":
	taps = []
	for i in range(1,len(sys.argv)):
		taps.append(int(sys.argv[i]))
	print('Taps: ', taps)
	order = taps[0]
	taps.pop(0)				# remove the first (highest order) tap, the way we do LFSR below you don't XOR that one
	print('Order: ', order)
	fills = generateFills(order)

	for fill in fills:
		print('')
		print('Fill:  ',fill)
		print('Output: ',generic(order, taps, fill))
		#print('Output: ',four(fill))
		exit(0)
