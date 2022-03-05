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

	// inputDumpFilename = getInputFileName_direct (".lammpstrj");
	// inputDataFilename = getInputFileName_direct (".data");

	if (argc == 3)
	{
		inputDumpFilename = argv[1];
		inputDataFilename = argv[2];
	}
	else
	{
		printf("\nREQUIRED ARGS:\n~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = dump filename\nargv[2] = data filename\n\n");
		exit (1);
	}

	inputConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputConfigFilename, 50, "bondRDF.config");

	inputFreevolumeConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputFreevolumeConfigFilename, 50, "freeVolume.config");

	inputVDWConfigFilename = (char *) malloc (50 * sizeof (char));
	snprintf (inputVDWConfigFilename, 50, "vdwsize.config");

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

	computeHBondCorrelation (inputDumpFile, nThreads);

	computeACFOfBondRDF (inputDumpFile);

	return 0;
}