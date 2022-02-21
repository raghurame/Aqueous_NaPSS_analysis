import csv
from time import sleep
import os
import numpy as n
from matplotlib import pyplot as plt
import statistics as s
from collections import Counter as c
import math as m

def showPlot (y, title, xLabel, yLabel):
	x = n.arange (1, len (y) + 1)
	plt.title (title)
	plt.xlabel (xLabel)
	plt.ylabel (yLabel)
	plt.plot (x, y)
	plt.show ()

def getOnlyPositive (inputData):
	for x in range (len (inputData)):
		if inputData[x] < 0:
			inputData[x] = 0

	return inputData

def getOnlyNegative (inputData):
	for x in range (len (inputData)):
		if inputData[x] > 0:
			inputData[x] = 0

	return inputData

def computeTotalArea (inputData):
	inputNumpy = n.array (inputData)
	# showPlot (inputNumpy, "Unmodified ACF plot", "Time lag", "ACF")
	totalArea = n.trapz (inputNumpy, dx = 1)
	return totalArea

def computePositiveArea (inputData):
	inputData = getOnlyPositive (inputData)
	inputNumpy = n.array (inputData)
	# showPlot (inputNumpy, "Positive area", "Time lag", "ACF")
	positiveArea = n.trapz (inputNumpy, dx = 1)
	return positiveArea

def computeNegativeArea (inputData):
	inputData = getOnlyNegative (inputData)
	inputNumpy = n.array (inputData)
	# showPlot (inputNumpy, "Negative area", "Time lag", "ACF")
	negativeArea = n.trapz (inputNumpy, dx = 1)
	return negativeArea

def processRDFfiles (inputFilename):
	column1 = []
	column2 = []

	with open ("bondRDF_processed/" + inputFilename) as inputRDFfile:
		for line in inputRDFfile:
			column1.append (float ((line.split (",")[0]).replace ("\n", "")))
			column2.append (float ((line.split (",")[1]).replace ("\n", "")))

	totalArea = computeTotalArea (column2)
	positiveArea = computePositiveArea (column2)
	negativeArea = computeNegativeArea (column2)

	with open ("RDF_AUC.output", "a") as outputFile:
		outputFile.write ("%.3f, %.3f, %.3f, %.3f\n" % (column1[0], positiveArea, negativeArea, totalArea))

	return column1[0], positiveArea, negativeArea, totalArea

def getDistribution (list_col1, list_posArea, list_negArea, list_totArea):
	nElements = len (list_totArea)
	maxPosValue = max (list_posArea)
	maxNegValue = max (list_negArea)
	maxTotValue = max (list_totArea)

	nBinsPos = 20
	nBinsNeg = 20
	nBinsTot = 20

	binSizePos = m.ceil (maxPosValue / nBinsPos)
	binSizeNeg = m.ceil (maxNegValue / nBinsNeg)
	binSizeTot = m.ceil (maxTotValue / nBinsTot)

	distTot_peak1 = [0] * nBinsTot
	distTot_peak2 = [0] * nBinsTot
	distTot_peak3 = [0] * nBinsTot
	distNeg_peak1 = [0] * nBinsNeg
	distNeg_peak2 = [0] * nBinsNeg
	distNeg_peak3 = [0] * nBinsNeg
	distPos_peak1 = [0] * nBinsPos
	distPos_peak2 = [0] * nBinsPos
	distPos_peak3 = [0] * nBinsPos

	# Distribution of total area
	binLow = 0
	for x in range (nBinsTot):
		binHigh = binLow + binSizeTot
		for y in range (len (list_totArea)):
			if list_totArea[y] < binHigh and list_totArea[y] > binLow:
				if list_col1[y] > 0 and list_col1[y] <= 1.9:
					distTot_peak1[x] += 1
				elif list_col1[y] > 1.9 and list_col1[y] <= 4:
					distTot_peak2[x] += 1
				elif list_col1[y] > 4 and list_col1[y] <= 6.6:
					distTot_peak3[x] += 1
		binLow = binHigh

	# Distribution of positive area
	binLow = 0
	for x in range (nBinsPos):
		binHigh = binLow + binSizePos
		for y in range (len (list_posArea)):
			if list_posArea[y] < binHigh and list_posArea[y] > binLow:
				if list_col1[y] > 0 and list_col1[y] <= 1.9:
					distPos_peak1[x] += 1
				elif list_col1[y] > 1.9 and list_col1[y] <= 4:
					distPos_peak2[x] += 1
				elif list_col1[y] > 4 and list_col1[y] <= 6.6:
					distPos_peak3[x] += 1
		binLow = binHigh

	for freq in distPos_peak1:
		print (freq)

	# Distribution of negative area
	binLow = 0
	for x in range (nBinsNeg):
		binHigh = binLow + binSizeNeg
		for y in range (len (list_negArea)):
			if list_negArea[y] < binHigh and list_negArea[y] > binLow:
				if list_col1[y] > 0 and list_col1[y] <= 1.9:
					distNeg_peak1[x] += 1
				elif list_col1[y] > 1.9 and list_col1[y] <= 4:
					distNeg_peak2[x] += 1
				elif list_col1[y] > 4 and list_col1[y] <= 6.6:
					distNeg_peak3[x] += 1
		binLow = binHigh

	distTotP1Numpy = n.array (distTot_peak1)
	distTotP2Numpy = n.array (distTot_peak2)
	distTotP3Numpy = n.array (distTot_peak3)
	distPosP1Numpy = n.array (distPos_peak1)
	distPosP2Numpy = n.array (distPos_peak2)
	distPosP3Numpy = n.array (distPos_peak3)
	distNegP1Numpy = n.array (distNeg_peak1)
	distNegP2Numpy = n.array (distNeg_peak2)
	distNegP3Numpy = n.array (distNeg_peak3)
	showPlot (distTotP1Numpy, "Distribution: Total area - Peak 1", "Area", "Freq")
	showPlot (distTotP2Numpy, "Distribution: Total area - Peak 2", "Area", "Freq")
	showPlot (distTotP3Numpy, "Distribution: Total area - Peak 3", "Area", "Freq")
	showPlot (distPosP1Numpy, "Distribution: Positive area - Peak 1", "Area", "Freq")
	showPlot (distPosP2Numpy, "Distribution: Positive area - Peak 2", "Area", "Freq")
	showPlot (distPosP3Numpy, "Distribution: Positive area - Peak 3", "Area", "Freq")
	showPlot (distNegP1Numpy, "Distribution: Negative area - Peak 1", "Area", "Freq")
	showPlot (distNegP2Numpy, "Distribution: Negative area - Peak 2", "Area", "Freq")
	showPlot (distNegP3Numpy, "Distribution: Negative area - Peak 3", "Area", "Freq")

def main ():
	list_col1 = []
	list_posArea = []
	list_negArea = []
	list_totArea = []
	header = []
	for dirs, dirname, files in os.walk ("bondRDF_processed/"):
		for file in files:
			col1, posArea, negArea, totArea = processRDFfiles (file)
			list_col1.append (col1)
			list_posArea.append (posArea)
			list_negArea.append (negArea)
			list_totArea.append (totArea)

	getDistribution (list_col1, list_posArea, list_negArea, list_totArea)

if __name__ == '__main__':
	os.system ("touch RDF_AUC.output")
	main()