import os
from time import sleep
import sys
import subprocess
import decimal

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def getNlines (dirs, filename):
	num_lines = sum (1 for line in open(dirs + "/" + filename))
	return num_lines

def computeACF (dirs, filename, nLines):
	print (nLines)
	with open (dirs + "/" + filename, "r") as inputFile:
		for line in inputFile:
			lineData = list (extract_numbers (line))
			print (lineData)
			sleep (1)

def main ():
	for dirs, dirname, files in os.walk ("."):
		for file in files:
			if ".oop" in file:
				nLines = getNlines (dirs, file)
				computeACF (dirs, file, nLines)

if __name__ == '__main__':
	main ()