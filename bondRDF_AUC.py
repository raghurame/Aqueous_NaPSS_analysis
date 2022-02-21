import csv
from time import sleep
import os
import numpy as n
from matplotlib import pyplot as plt

def showPlot (y, title, xLabel, yLabel):
	x = n.arange (1, len (y) + 1)
	plt.title (title)
	plt.xlabel (xLabel)
	plt.ylabel (yLabel)
	plt.plot (x, y)
	plt.show ()

def getOnlyPositive (inputData):
	for item in inputData:
		if item < 0:
			item = 0

	return inputData

def getOnlyNegative (inputData):
	for item in inputData:
		if item > 0:
			item = 0

	return inputData

def computeTotalArea (inputData):
	inputNumpy = n.array (inputData)
	showPlot (inputNumpy, "Unmodified ACF plot", "Time lag", "ACF")
	totalArea = n.trapz (inputNumpy, dx = 1)
	return totalArea

def computePositiveArea (inputData):
	inputData = getOnlyPositive (inputData)
	inputNumpy = n.array (inputData)
	showPlot (inputNumpy, "Positive area", "Time lag", "ACF")
	positiveArea = n.trapz (inputNumpy, dx = 1)
	return positiveArea

def computeNegativeArea (inputData):
	inputData = getOnlyNegative (inputData)
	inputNumpy = n.array (inputData)
	showPlot (inputNumpy, "Negative area", "Time lag", "ACF")
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

	print ("Positive area: {}\nNegative area: {}\nTotal area: {}\n".format (positiveArea, negativeArea, totalArea))

def main ():
	header = []
	for dirs, dirname, files in os.walk ("bondRDF_processed/"):
		for file in files:
			processRDFfiles (file)

if __name__ == '__main__':
	main()