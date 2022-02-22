import csv
from time import sleep
import os
import numpy as n
from matplotlib import pyplot as plt
import statistics as s
from collections import Counter as c
import math as m

def storeFileValues (targetDirectory, filename, inputData):
	with open (targetDirectory + filename, "r") as inputFile:
		for line in inputFile:
			a, b, c, d, e = line.split (" ")
			inputData[a] = [b, c, d, e]

	return inputData

def findMSD (msd, inputData):
	

def processPeak (targetDirectory, fileNameTemplate):
	inputData = {}
	msd = []
	for dirs, dirname, files in os.walk (targetDirectory):
		for file in files:
			if (fileNameTemplate in file):
				print (file)
				inputData = storeFileValues (targetDirectory, file, inputData)
				msd = findMSD (msd, inputData)
				print (inputData)
				print (msd)
				exit (1)

def main ():
	processPeak ("MSD_logs/", "_1.msd")
	processPeak ("MSD_logs/", "_2.msd")
	processPeak ("MSD_logs/", "_3.msd")

if __name__ == '__main__':
	main()