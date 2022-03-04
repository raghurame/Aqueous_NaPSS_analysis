import os
from time import sleep
import sys
import subprocess
import decimal
from pathlib import Path
from tqdm import tqdm

'''
Arguments to pass:
~~~~~~~~~~~~~~~~~

argv[0] = ./program
argv[1] = main LAMMPS traj dump filename
argv[2] = threshold distance to consider from the first timeframe

'''

def extract_numbers (line):
	lineString = line.replace ("\t", ",").replace (" ", ",").replace ("\n", "")
	for item in lineString.split (","):
		try:
			yield decimal.Decimal (item)
		except:
			pass

def readLogFile (logFilename, mainData):
	thresholdDistance = float (sys.argv[2])
	lineNumber = 1
	firstTimeframe = 0
	lineData = []

	if (len (mainData) == 0):
		firstTimeframe = 1

	with open (logFilename, "r") as inputLogFile:
		for line in tqdm (inputLogFile, desc = "reading {}".format (logFilename)):
			lineData = list (extract_numbers (line))
			if (firstTimeframe == 1):
				mainData[lineNumber] = []
				
			mainData[lineNumber].append (lineData[4])
			lineNumber += 1

	return mainData

def main ():
	mainDumpfile = sys.argv[1]
	isTimestep = 0
	mainData = {}

	with open (mainDumpfile, "r") as inputFile:
		for line in inputFile:
			if (isTimestep == 1):
				line = line.replace ("\n", "")
				logFilename = ("bondRDF_logs/" + line + ".rdf")
				path = Path (logFilename)
				if (path.is_file ()):
					mainData = readLogFile (logFilename, mainData)
				else:
					print (logFilename + " not found")
				isTimestep = 0

			if ("ITEM: TIMESTEP" in line):
				isTimestep = 1

	for key in mainData:
		print (key, mainData[key])
		sleep (1)

if __name__ == '__main__':
	main ()