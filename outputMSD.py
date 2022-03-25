import csv
from time import sleep
import os
import numpy as n
from matplotlib import pyplot as plt
import statistics as s
from collections import Counter as c
import math as m
from tqdm import tqdm
import sys

def storeFileValues (targetDirectory, filename):
	inputData = []
	with open (targetDirectory + filename, "r") as inputFile:
		for line in inputFile:
			line.replace ("\n", "")
			a, b, c, d, e = line.split (" ")
			# inputData array structure:
			# currentDumpstep, distance between SO and OH, X, Y, Z.
			inputData.append ([int (a), float (b), float (c), float (d), float (e)])

	return inputData

def findMSD (msd, inputData):
	# msd dict structure:
	# delT, [sum of displacementSquare, number of elements in sum]
	for t0 in inputData:
		for tx in inputData:
			if ((t0[0] != tx[0]) and (tx[0] > t0[0])):
				delT = tx[0] - t0[0]
				displacement = m.dist ([t0[2], t0[3], t0[4]], [tx[2], tx[3], tx[4]])
				displacementSquare = displacement**2

				try:
					msd[delT][0] += displacementSquare
					msd[delT][1] += 1
				except:
					msd[delT] = [displacementSquare, 1]
			
	# for key in sorted (msd.keys ()):
	# 	print (key, msd[key][0], float (msd[key][0] / msd[key][1]), msd[key][1])

	return msd

def processPeak (targetDirectory, fileNameTemplate, delT_dump):
	inputData = []
	msd = {}
	for dirs, dirname, files in os.walk (targetDirectory):
		for file in tqdm (files, desc = "reading peak: {}".format (fileNameTemplate.replace ("_", "").replace (".msd", ""))):
			if (fileNameTemplate in file):
				inputData = storeFileValues (targetDirectory, file)
				msd = findMSD (msd, inputData)

	logX = []
	logY = []

	with open ("msd_log" + fileNameTemplate, "w") as outputMSD:
		for key in sorted (msd.keys ()):
			outputMSD.write ("{}, {}\n".format (m.log (key * int (delT_dump), 10), m.log (float (msd[key][0] / msd[key][1]), 10)))
			logX.append (m.log (key * int (delT_dump), 10))
			logY.append (m.log (float (msd[key][0] / msd[key][1]), 10))

	nplogX = n.array (logX)
	nplogY = n.array (logY)

	slope_intercept = n.polyfit (nplogX, nplogY, 1)
	print (slope_intercept)

	with open ("diffusionConstant_log" + fileNameTemplate, "w") as outputDiffusionConstant:
		outputDiffusionConstant.write (str (slope_intercept[0]/6) + "\n")

	logX = []
	logY = []

	with open ("msd" + fileNameTemplate, "w") as outputMSD:
		for key in sorted (msd.keys ()):
			outputMSD.write ("{}, {}\n".format (key * int (delT_dump), float (msd[key][0] / msd[key][1])))
			logX.append (float (key * int (delT_dump)))
			logY.append (float (msd[key][0] / msd[key][1]))

	lenX = m.ceil (len (logX) / 2)
	lenY = m.ceil (len (logY) / 2)

	for x in range (0, (lenX * 2)):
		print (logX[x], logY[x])
		sleep (1)

	nplogX = n.array (logX[:lenX])
	nplogY = n.array (logY[:lenY])

	slope_intercept = n.polyfit (nplogX, nplogY, 1)
	print (slope_intercept)

	with open ("diffusionConstant" + fileNameTemplate, "w") as outputDiffusionConstant:
		outputDiffusionConstant.write (str (slope_intercept[0]/6) + "\n")

def main (delT_dump):
	processPeak ("MSD_logs/", "_1.msd", delT_dump)
	processPeak ("MSD_logs/", "_2.msd", delT_dump)
	processPeak ("MSD_logs/", "_3.msd", delT_dump)
	processPeak ("MSD_logs/", "_4.msd", delT_dump)

if __name__ == '__main__':
	main(sys.argv[1])