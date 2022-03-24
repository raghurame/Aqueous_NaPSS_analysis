from time import sleep
import sys
import statistics as s

def main (inputFilename):
	NfirstToSecond = []
	NsecondToThird = []
	NthirdToFourth = []
	NMolsZeroth = []
	NMolsFirst = []
	NMolsSecond = []
	NMolsThird = []
	NMolsFourth = []
	lineData = []
	first_fraction = []
	second_fraction = []
	third_fraction = []

	with open (inputFilename, "r") as inputfile:
		for line in inputfile:
			if (line[0] != '#'):
				lineData = line.split (",")

				NfirstToSecond.append (int (lineData[0]))
				NsecondToThird.append (int (lineData[1]))
				NthirdToFourth.append (int (lineData[2]))

				NMolsZeroth.append (int (lineData[3]))
				NMolsFirst.append (int (lineData[4]))
				NMolsSecond.append (int (lineData[5]))
				NMolsThird.append (int (lineData[6]))
				NMolsFourth.append (int (lineData[7]))

				first_fraction.append (float (lineData[0]) / float (lineData[4]))
				second_fraction.append (float (lineData[1]) / float (lineData[5]))
				third_fraction.append (float (lineData[2]) / float (lineData[6]))

	with open ("hBondNetwork.stats", "w") as outputFile:
		outputFile.write ("hBond density in first layer: {} (+/- {})\nhBond density in second layer: {} (+/- {})\nhBond density in third layer: {} (+/- {})\n".format (s.mean (first_fraction), s.stdev (first_fraction), s.mean (second_fraction), s.stdev (second_fraction), s.mean (third_fraction), s.stdev (third_fraction)))

	print ("hBond density in first layer: {} (+/- {})\nhBond density in second layer: {} (+/- {})\nhBond density in third layer: {} (+/- {})\n".format (s.mean (first_fraction), s.stdev (first_fraction), s.mean (second_fraction), s.stdev (second_fraction), s.mean (third_fraction), s.stdev (third_fraction)))

if __name__ == '__main__':
	if (len (sys.argv) != 2):
		print ("REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\nargv[0] = program\nargv[1] = input filename\n")
		exit (1)
	main(sys.argv[1])