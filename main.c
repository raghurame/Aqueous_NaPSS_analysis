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
	
	// nThreads = 1; // Comment this line if enough processors are available 

	FILE *inputDumpFile, *inputDataFile, *inputConfigFile, *inputFreevolumeConfigFile, *inputVDWConfigFile, *inputHBondConfigFile;
	char *inputDumpFilename, *inputDataFilename, *inputConfigFilename, *inputFreevolumeConfigFilename, *inputVDWConfigFilename, *inputHBondConfigFilename;

	// Uncomment the following two lines if only one LAMMPSTRJ and DATA files are available
	// inputDumpFilename = getInputFileName_direct (".lammpstrj");
	// inputDataFilename = getInputFileName_direct (".data");

	// If the above two lines are uncommented, then the following if-else condition must be commented
	if (argc == 4)
	{
		inputDumpFilename = argv[1];
		inputDataFilename = argv[2];
		nThreads = atoi (argv[3]);
	}
	else
	{
		printf("\nREQUIRED ARGS:\n~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = dump filename\nargv[2] = data filename\nargv[3] = no. of processors\n\n");
		exit (1);
	}

	int bondRDF_config_status, freeVolume_config_status, vdwsize_config_status, hBond_config_status, msd_config_status, hBond_threshold_status, freeVolumeVars_config_status;

	printf("\nChecking necessary config files...\n\n");
	bondRDF_config_status = verifyConfigFiles ("bondRDF.config");
	freeVolume_config_status = verifyConfigFiles ("freeVolume.config");
	vdwsize_config_status = verifyConfigFiles ("vdwsize.config");
	hBond_config_status = verifyConfigFiles ("hBond.config");
	msd_config_status = verifyConfigFiles ("msd.config");
	hBond_threshold_status = verifyConfigFiles ("hBond.threshold");
	freeVolumeVars_config_status = verifyConfigFiles ("freeVolumeVars.config");

	if ((bondRDF_config_status + freeVolume_config_status + vdwsize_config_status + hBond_config_status + msd_config_status + hBond_threshold_status + freeVolumeVars_config_status) > 0)
	{
		printf("\nVital files missing. Program will exit now !\n\n");
		exit (1);
	}

	inputDumpFile = fopen (inputDumpFilename, "r");
	inputDataFile = fopen (inputDataFilename, "r");
	inputConfigFile = fopen ("bondRDF.config", "r");
	inputFreevolumeConfigFile = fopen ("freeVolume.config", "r");
	inputVDWConfigFile = fopen ("vdwsize.config", "r");
	inputHBondConfigFile = fopen ("hBond.config", "r");

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

	free (atoms);
	free (bonds);
	free (angles);
	free (dihedrals);
	free (impropers);

	free (inputVectors);
	free (freeVolumeconfig);
	free (vwdSize);
	free (HBondAtoms);

	return 0;
}