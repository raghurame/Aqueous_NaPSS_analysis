#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <omp.h>
#include "structDefinitions.h"
#include "aqueousnapss.h"
#include "fileHandling.h"
#include "readInputFile.h"
#include "hBondCorrelation.h"
#include "bondRDF_ACF.h"

int main (int argc, char const *argv[])
{
	system ("mkdir logs");
	system ("mkdir bondRDF_logs");
	system ("mkdir hBondNetwork_logs");

	long number_of_processors = sysconf(_SC_NPROCESSORS_ONLN);
	int nThreads = (int) number_of_processors - 1;

	FILE *inputDumpFile, *inputDataFile, *inputConfigFile, *inputFreevolumeConfigFile, *inputVDWConfigFile, *inputHBondConfigFile;
	char *inputDumpFilename, *inputDataFilename, *inputConfigFilename, *inputFreevolumeConfigFilename, *inputVDWConfigFilename, *inputHBondConfigFilename;

	// printf("%s\n", "Looking for LAMMPS trajectory file...");
	// inputDumpFilename = getInputFileName ();
	inputDumpFilename = argv[1];
	// printf("%s\n", "Looking for LAMMPS data file...");
	// inputDataFilename = getInputFileName ();
	inputDataFilename = argv[2];
	// printf("%s\n", "Looking for input config file (for OOP/bondRDF calculations)...");
	// inputConfigFilename = getInputFileName ();
	inputConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputConfigFilename, 50, "bondRDF.config");
	// printf("%s\n", "Looking for config file for free volume calculations...");
	// inputFreevolumeConfigFilename = getInputFileName ();
	inputFreevolumeConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputFreevolumeConfigFilename, 50, "freeVolume.config");
	// printf("%s\n", "Looking for config file containing VWD radii...");
	// inputVDWConfigFilename = getInputFileName ();
	inputVDWConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputVDWConfigFilename, 50, "vdwsize.config");
	// printf("%s\n", "Looking for config file to compute H bonding...");
	// inputHBondConfigFilename = getInputFileName ();
	inputHBondConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputHBondConfigFilename, 50, "hBond.config");

	inputDumpFile = fopen (inputDumpFilename, "r");
	inputDataFile = fopen (inputDataFilename, "r");
	inputConfigFile = fopen (inputConfigFilename, "r");
	inputFreevolumeConfigFile = fopen (inputFreevolumeConfigFilename, "r");
	inputVDWConfigFile = fopen (inputVDWConfigFilename, "r");
	inputHBondConfigFile = fopen (inputHBondConfigFilename, "r");

	DATA_ATOMS *atoms;
	DATA_BONDS *bonds;
	DATA_ANGLES *angles;
	DATA_DIHEDRALS *dihedrals;
	DATA_IMPROPERS *impropers;

	DATAFILE_INFO datafile;
	datafile = readData (inputDataFile, &atoms, &bonds, &angles, &dihedrals, &impropers);

	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);

	CONFIG *inputVectors, *freeVolumeconfig, *vwdSize, *HBondAtoms;
	NLINES_CONFIG entries;
	int nLines_inputVectors, nLines_freeVolumeconfig, nLines_vwdSize, nLines_HBondAtoms;

	inputVectors = readConfig (inputConfigFile, &nLines_inputVectors);
	freeVolumeconfig = readConfig (inputFreevolumeConfigFile, &nLines_freeVolumeconfig);
	vwdSize = readVWDRadius (inputVDWConfigFile, &nLines_vwdSize);
	HBondAtoms = readConfig (inputHBondConfigFile, &nLines_HBondAtoms);

	entries.nLines_inputVectors = nLines_inputVectors; entries.nLines_freeVolumeconfig = nLines_freeVolumeconfig; entries.nLines_vwdSize = nLines_vwdSize; entries.nLines_HBondAtoms = nLines_HBondAtoms;

	processLAMMPSTraj (inputDumpFile, datafile, bonds, inputVectors, freeVolumeconfig, vwdSize, entries, nThreads);

	// Checking h-bond lifetime correlation function from the saved logfiles
	computeHBondCorrelation (inputDumpFile, nThreads);

	// Computing autocorrelation of bondRDF distances
	computeACFOfBondRDF (inputDumpFile);

	return 0;
}