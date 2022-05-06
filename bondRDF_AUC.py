import csv
from time import sleep
import os
import numpy as n
from matplotlib import pyplot as plt
import statistics as s
from collections import Counter as c
import math as m
from tqdm import tqdm

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
	totalArea = n.trapz (inputNumpy, dx = 0.5)
	return totalArea

def computePositiveArea (inputData):
	inputData = getOnlyPositive (inputData)
	inputNumpy = n.array (inputData)
	# showPlot (inputNumpy, "Positive area", "Time lag", "ACF")
	positiveArea = n.trapz (inputNumpy, dx = 0.5)
	return positiveArea

def computeNegativeArea (inputData):
	inputData = getOnlyNegative (inputData)
	inputNumpy = n.array (inputData)
	# showPlot (inputNumpy, "Negative area", "Time lag", "ACF")
	negativeArea = n.trapz (inputNumpy, dx = 0.5)
	return negativeArea

def keepFirstPositive (column2):
	column2_firstPositive = []

	for x in range (0, len (column2), 1):
		if (column2[x] <= 0):
			for y in range (0, x, 1):
				column2_firstPositive.append (column2[y])
			return column2_firstPositive

def processRDFfiles (inputFilename):
	column1 = []
	column2 = []

	with open ("bondRDF_processed/" + inputFilename) as inputRDFfile:
		for line in inputRDFfile:
			column1.append (float ((line.split (",")[0]).replace ("\n", "")))
			column2.append (float ((line.split (",")[1]).replace ("\n", "")))

	column2_firstPositive = keepFirstPositive (column2)

	totalArea = computeTotalArea (column2)
	positiveArea = computePositiveArea (column2_firstPositive)
	negativeArea = computeNegativeArea (column2)

	with open ("RDF_AUC.output", "a") as outputFile:
		outputFile.write ("%.3f, %.3f, %.3f, %.3f\n" % (column1[0], positiveArea, negativeArea, totalArea))

	return column1[0], positiveArea, negativeArea, totalArea

def normalizeDistribution (rawDistribution):
	maxValue = max (rawDistribution)

	if (maxValue == 0):
		return rawDistribution
	else:
		for x in range (len (rawDistribution)):
			rawDistribution[x] = rawDistribution[x] / maxValue

	return rawDistribution

def multiPlot (list1, label1, list2, label2, list3, label3, list4, label4, xLabel, yLabel, title):
	x = n.arange (1, len (list1) + 1)

	plt.plot (x, list1, label = label1)
	plt.plot (x, list2, label = label2)
	plt.plot (x, list3, label = label3)
	plt.plot (x, list4, label = label4)

	plt.xlabel (xLabel)
	plt.ylabel (yLabel)
	plt.title (title)

	plt.legend ()

	plt.show ()

def getDistribution (list_col1, list_posArea, list_negArea, list_totArea):

	peak1_min = 2.5
	peak1_max = 5.0
	peak2_min = 5.0
	peak2_max = 8.6
	peak3_min = 8.6
	peak3_max = 12.2
	peak4_min = 12.2

	nElements = len (list_totArea)
	maxPosValue = max (list_posArea)
	maxNegValue = max (list_negArea)
	maxTotValue = max (list_totArea)

	nBinsPos = 50
	nBinsNeg = 50
	nBinsTot = 50

	binSizePos = m.ceil (maxPosValue / nBinsPos)
	binSizeNeg = m.ceil (maxNegValue / nBinsNeg)
	binSizeTot = m.ceil (maxTotValue / nBinsTot)

	distTot_peak1 = [0] * nBinsTot
	distTot_peak2 = [0] * nBinsTot
	distTot_peak3 = [0] * nBinsTot
	distTot_peak4 = [0] * nBinsTot
	distNeg_peak1 = [0] * nBinsNeg
	distNeg_peak2 = [0] * nBinsNeg
	distNeg_peak3 = [0] * nBinsNeg
	distNeg_peak4 = [0] * nBinsNeg
	distPos_peak1 = [0] * nBinsPos
	distPos_peak2 = [0] * nBinsPos
	distPos_peak3 = [0] * nBinsPos
	distPos_peak4 = [0] * nBinsPos

	# Distribution of total area
	binLow = 0
	for x in tqdm (range (nBinsTot), desc = "Distribution of total area"):
		binHigh = binLow + binSizeTot
		for y in range (len (list_totArea)):
			if list_totArea[y] < binHigh and list_totArea[y] > binLow:
				if list_col1[y] > peak1_min and list_col1[y] <= peak1_max:
					distTot_peak1[x] += 1
				elif list_col1[y] > peak2_min and list_col1[y] <= peak2_max:
					distTot_peak2[x] += 1
				elif list_col1[y] > peak3_min and list_col1[y] <= peak3_max:
					distTot_peak3[x] += 1
				elif list_col1[y] > peak4_min:
					distTot_peak4[x] += 1
		binLow = binHigh

	# Distribution of positive area
	binLow = 0
	for x in tqdm (range (nBinsPos), desc = "Distribution of positive area"):
		binHigh = binLow + binSizePos
		for y in range (len (list_posArea)):
			if list_posArea[y] < binHigh and list_posArea[y] > binLow:
				if list_col1[y] > peak1_min and list_col1[y] <= peak1_max:
					distPos_peak1[x] += 1
				elif list_col1[y] > peak2_min and list_col1[y] <= peak2_max:
					distPos_peak2[x] += 1
				elif list_col1[y] > peak3_min and list_col1[y] <= peak3_max:
					distPos_peak3[x] += 1
				elif list_col1[y] > peak4_min:
					distPos_peak4[x] += 1
		binLow = binHigh

	# Distribution of negative area
	binLow = 0
	for x in tqdm (range (nBinsNeg), desc = "Distribution of negative area"):
		binHigh = binLow + binSizeNeg
		for y in range (len (list_negArea)):
			if list_negArea[y] < binHigh and list_negArea[y] > binLow:
				if list_col1[y] > peak1_min and list_col1[y] <= peak1_max:
					distNeg_peak1[x] += 1
				elif list_col1[y] > peak2_min and list_col1[y] <= peak2_max:
					distNeg_peak2[x] += 1
				elif list_col1[y] > peak3_min and list_col1[y] <= peak3_max:
					distNeg_peak3[x] += 1
				elif list_col1[y] > peak4_min:
					distNeg_peak4[x] += 1
		binLow = binHigh

	# distTot_peak1 = normalizeDistribution (distTot_peak1)
	# distTot_peak2 = normalizeDistribution (distTot_peak2)
	# distTot_peak3 = normalizeDistribution (distTot_peak3)
	# distTot_peak4 = normalizeDistribution (distTot_peak4)

	distPos_peak1 = normalizeDistribution (distPos_peak1)
	distPos_peak2 = normalizeDistribution (distPos_peak2)
	distPos_peak3 = normalizeDistribution (distPos_peak3)
	distPos_peak4 = normalizeDistribution (distPos_peak4)

	with open ("bondRDF_AUC.output", "w") as bondRDFAUCfile:
		for x in range (0, nBinsPos, 1):
			bondRDFAUCfile.write ("{}, {}, {}, {}\n".format (distPos_peak1[x], distPos_peak2[x], distPos_peak3[x], distPos_peak4[x]))

	# distNeg_peak1 = normalizeDistribution (distNeg_peak1)
	# distNeg_peak2 = normalizeDistribution (distNeg_peak2)
	# distNeg_peak3 = normalizeDistribution (distNeg_peak3)
	# distNeg_peak4 = normalizeDistribution (distNeg_peak4)

	# distTotP1Numpy = n.array (distTot_peak1)
	# distTotP2Numpy = n.array (distTot_peak2)
	# distTotP3Numpy = n.array (distTot_peak3)
	# distTotP4Numpy = n.array (distTot_peak4)

	distPosP1Numpy = n.array (distPos_peak1)
	distPosP2Numpy = n.array (distPos_peak2)
	distPosP3Numpy = n.array (distPos_peak3)
	distPosP4Numpy = n.array (distPos_peak4)

	# distNegP1Numpy = n.array (distNeg_peak1)
	# distNegP2Numpy = n.array (distNeg_peak2)
	# distNegP3Numpy = n.array (distNeg_peak3)
	# distNegP4Numpy = n.array (distNeg_peak4)

	# showPlot (distTotP1Numpy, "Distribution: Total area - Peak 1", "Area", "Freq")
	# showPlot (distTotP2Numpy, "Distribution: Total area - Peak 2", "Area", "Freq")
	# showPlot (distTotP3Numpy, "Distribution: Total area - Peak 3", "Area", "Freq")

	# showPlot (distPosP1Numpy, "Distribution: Positive area - Peak 1", "Area", "Freq")
	# showPlot (distPosP2Numpy, "Distribution: Positive area - Peak 2", "Area", "Freq")
	# showPlot (distPosP3Numpy, "Distribution: Positive area - Peak 3", "Area", "Freq")
	# showPlot (distPosP4Numpy, "Distribution: Positive area - Bulk", "Area", "Freq")

	# showPlot (distNegP1Numpy, "Distribution: Negative area - Peak 1", "Area", "Freq")
	# showPlot (distNegP2Numpy, "Distribution: Negative area - Peak 2", "Area", "Freq")
	# showPlot (distNegP3Numpy, "Distribution: Negative area - Peak 3", "Area", "Freq")

	multiPlot (distPosP1Numpy, "peak 1", distPosP2Numpy, "peak 2", distPosP3Numpy, "peak 3", distPosP4Numpy, "peak 4", "Area Distribution", "Freq", "ACF of Dipole Distances - Positive Area")


def main ():
	list_col1 = []
	list_posArea = []
	list_negArea = []
	list_totArea = []
	header = []
	for dirs, dirname, files in os.walk ("bondRDF_processed/"):
		for file in tqdm (files, desc = "reading bondRDF files"):
			col1, posArea, negArea, totArea = processRDFfiles (file)
			list_col1.append (col1)
			list_posArea.append (posArea)
			list_negArea.append (negArea)
			list_totArea.append (totArea)

	getDistribution (list_col1, list_posArea, list_negArea, list_totArea)

if __name__ == '__main__':
	os.system ("touch RDF_AUC.output")
	main()