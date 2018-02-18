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

import numpy as np
import math
import sys


def ASCIIToDecimal(ascii):
	decimal = np.zeros(len(ascii), dtype = np.int)
	for index, char in enumerate(ascii):
		decimal[index] = ord(char)
	return decimal


# This reads the binary from MSB left to LSB right (big endian)
# This will create N number of ASCII characters = binary length / 8 bits
# If there are 'leftover' bits, b/c the length isn't evenly divisable by 8, 
# then the extra bits on the right side will be thrown out.
def BinaryToASCII(binary):
	if ((len(binary)%8) != 0):
		print("Warning: number of bits does not evenly divide into 8 bit characters. Will read from left and throw out leftover.")
		# this means there are leftover bits, maybe error instead
	
	numchar = int(math.floor(len(binary)/8))

	ascii = []
	j = 0
	for i in range(numchar):
		decimal = binary[j]*128 + binary[j+1]*64 + binary[j+2]*32 + binary[j+3]*16 + binary[j+4]*8 + binary[j+5]*4 + binary[j+6]*2 + binary[j+7]*1
		ascii.append(chr(decimal))
		j = j + 8
	
	return ascii
	
	
# This reads the binary from MSB left to LSB right (big endian)
# This will create N number of ASCII decimals = binary length / 8 bits
# If there are 'leftover' bits, b/c the length isn't evenly divisable by 8, 
# then the extra bits on the right side will be thrown out.
def BinaryToDecimalASCII(binary):
	if ((len(binary)%8) != 0):
		print("Warning: number of bits does not evenly divide into 8 bit characters. Will read from left and throw out leftover.")
		# this means there are leftover bits, maybe error instead
	
	numchar = int(math.floor(len(binary)/8))
	decimal_ascii = np.zeros((numchar), dtype = np.int)

	j = 0
	for i in range(numchar):
		decimal = binary[j]*128 + binary[j+1]*64 + binary[j+2]*32 + binary[j+3]*16 + binary[j+4]*8 + binary[j+5]*4 + binary[j+6]*2 + binary[j+7]*1
		decimal_ascii[i] = decimal
		j = j + 8
	
	return decimal_ascii

	
# This reads the binary from MSB left to LSB right (big endian)
# This will create N number of Hex characters = binary length / 4 bits
# If there are 'leftover' bits, b/c the length isn't evenly divisable by 4, 
# then the extra bits on the right side will be thrown out.
def BinaryToHex(binary):
	if ((len(binary)%4) != 0):
		print("Warning: number of bits does not evenly divide into 4 bit characters. Will read from left and throw out leftover.")
		# this means there are leftover bits, maybe error instead
	
	numchar = int(math.floor(len(binary)/4))
	
	hexnum = []
	j = 0
	for i in range(numchar):
		decimal = binary[j]*8 + binary[j+1]*4 + binary[j+2]*2 + binary[j+3]*1;
		hexnum.append(hex(decimal)[2:])
		j = j + 4
	
	return hexnum
	
	
def DecimalToASCII(decimal):
	ascii = []
	for num in decimal:
		ascii.append(chr(num))
	return ascii

	
# because we expect ASCII decimals [0:255], this creates 8 bit binary lists
def DecimalASCIIToBinary(decimal):
	binary = np.zeros((len(decimal)*8), dtype = np.int)
	j = 0
	for char in decimal:
		if (char > 255):
			#print("Warning: this function only accepts decimal values that fit within 8 bits, in the range 0:255. This value will be translated to zero")
			continue
		elif (char < 0):
			#print("Warning: this function does not accept negative numbers. This value will be translated to zero")
			continue
			
		np.put(binary, [j, j+1, j+2, j+3, j+4, j+5, j+6, j+7], list(map(int,format(char,'08b'))))
		j = j + 8
	
	return binary

	
def HexToBinary(hexnum):
	binary = np.zeros((len(hexnum)*4), dtype = np.int)
	j = 0
	for char in hexnum:
		np.put(binary, [j, j+1, j+2, j+3], list(map(int,format(int(char,16),'04b'))))
		j = j + 4
	
	return binary


def shiftLeftDecimalASCII(decimal, shift):
	shifted_decimal = np.zeros(len(decimal), dtype = np.int)

	for i in range(len(decimal)):
		# lower case shift
		if (96 < decimal[i] and decimal[i] < 123):
			shifted_decimal[i] = decimal[i] - shift
			if (shifted_decimal[i] < 97):
				shifted_decimal[i] = 123 - (97-shifted_decimal[i])
				
		# upper case shift
		elif (64 < decimal[i] and decimal[i] < 91):
			shifted_decimal[i] = decimal[i] - shift
			if (shifted_decimal[i] < 65):
				shifted_decimal[i] = 91 - (65-shifted_decimal[i])

		else:
			#print("Warning: this function only supports shifting alphabet characters A-Z (65:90) and a-z (97:122). This number will be shifted using the non alphabet ASCII map (ignoring A-Z a-z).")
			shifted_decimal[i] = decimal[i] - shift
			if (96 < shifted_decimal[i] and shifted_decimal[i] < 123):
				shifted_decimal[i] = 97 - (123 - shifted_decimal[i])
			if (64 < shifted_decimal[i] and shifted_decimal[i] < 91):
				shifted_decimal[i] = 65 - (91 - shifted_decimal[i])
			if (shifted_decimal[i] < 0):
				shifted_decimal[i] = (255 - shifted_decimal[i])

	return shifted_decimal
	
	
def shiftRightDecimalASCII(decimal, shift):
	shifted_decimal = np.zeros(len(decimal), dtype = np.int)

	for i in range(len(decimal)):
		# lower case shift
		if (96 < decimal[i] and decimal[i] < 123):
			shifted_decimal[i] = decimal[i] + shift
			if (shifted_decimal[i] > 122):
				shifted_decimal[i] = 96 + shifted_decimal[i]%122
				
		# upper case shift
		elif (64 < decimal[i] and decimal[i] < 91):
			shifted_decimal[i] = decimal[i] + shift
			if (shifted_decimal[i] > 90):
				shifted_decimal[i] = 64 + shifted_decimal[i]%90

		else:
			#print("Warning: this function only supports shifting alphabet characters A-Z (65:90) and a-z (97:122). This number will be shifted using the non alphabet ASCII map (ignoring A-Z a-z).")
			shifted_decimal[i] = decimal[i] + shift
			if (64 < shifted_decimal[i] and shifted_decimal[i] < 91):
				shifted_decimal[i] = 90 + (shifted_decimal[i] - 64)
			if (96 < shifted_decimal[i] and shifted_decimal[i] < 123):
				shifted_decimal[i] = 122 + (shifted_decimal[i] - 96)
			if (255 < shifted_decimal[i]):
				shifted_decimal[i] = (shifted_decimal[i] - 255)

	return shifted_decimal

	
if __name__ == "__main__":
	print('No default functionality')
