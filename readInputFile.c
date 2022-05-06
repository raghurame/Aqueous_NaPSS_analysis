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
#include "readInputFile.h"

DATAFILE_INFO readData (FILE *input, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers)
{
	printf("Reading LAMMPS data file...\n");
	rewind (input);

	int isAtomLine = 0, /*nAtoms = -1,*/ nAtomLine = 0;
	int isBondLine = 0, /*nBonds = -1,*/ nBondLine = 0;
	int isAngleLine = 0, /*nAngles = -1,*/ nAngleLine = 0;
	int isDihedralLine = 0, /*nDihedrals = -1,*/ nDihedralLine = 0;
	int isImproperLine = 0, /*nImpropers = -1,*/ nImproperLine = 0;
	int printHeaderInfo = 1;

	DATAFILE_INFO datafile;
	datafile.nAtoms = -1;
	datafile.nBonds = -1;
	datafile.nAngles = -1;
	datafile.nDihedrals = -1;
	datafile.nImpropers = -1;

	char lineString[1000];

	*atoms = NULL;
	*bonds = NULL;
	*angles = NULL;
	*dihedrals = NULL;
	*impropers = NULL;

	while ((fgets (lineString, 1000, input) != NULL))
	{
		if (strstr (lineString, "atoms"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtoms);
			fprintf(stdout, "nAtoms detected: %d\n", datafile.nAtoms);
			fflush (stdout);
			(*atoms) = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));
		}

		if (strstr (lineString, "bonds"))
		{
			sscanf (lineString, "%d \n", &datafile.nBonds);
			fprintf(stdout, "nBonds detected: %d\n", datafile.nBonds);
			fflush (stdout);
			(*bonds) = (DATA_BONDS *) malloc (datafile.nBonds * sizeof (DATA_BONDS));
		}

		if (strstr (lineString, "angles"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngles);
			fprintf(stdout, "nAngles detected: %d\n", datafile.nAngles);
			fflush (stdout);
			(*angles) = (DATA_ANGLES *) malloc (datafile.nAngles * sizeof (DATA_ANGLES));
		}

		if (strstr (lineString, "dihedrals"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedrals);
			fprintf(stdout, "nDihedrals detected: %d\n", datafile.nDihedrals);
			fflush (stdout);
			(*dihedrals) = (DATA_DIHEDRALS *) malloc (datafile.nDihedrals * sizeof (DATA_DIHEDRALS));
		}

		if (strstr (lineString, "impropers"))
		{
			sscanf (lineString, "%d \n", &datafile.nImpropers);
			fprintf(stdout, "nImpropers detected: %d\n", datafile.nImpropers);
			fflush (stdout);
			(*impropers) = (DATA_IMPROPERS *) malloc (datafile.nImpropers * sizeof (DATA_IMPROPERS));
		}

		if (strstr (lineString, "atom types"))
		{
			sscanf (lineString, "%d \n", &datafile.nAtomTypes);
			fprintf(stdout, "%s%d\n", "nAtomTypes detected: ", datafile.nAtomTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "bond types"))
		{
			sscanf (lineString, "%d \n", &datafile.nBondTypes);
			fprintf(stdout, "%s%d\n", "nBondTypes detected: ", datafile.nBondTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "angle types"))
		{
			sscanf (lineString, "%d \n", &datafile.nAngleTypes);
			fprintf(stdout, "%s%d\n", "nAngleTypes detected: ", datafile.nAngleTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "dihedral types"))
		{
			sscanf (lineString, "%d \n", &datafile.nDihedralTypes);
			fprintf(stdout, "%s%d\n", "nDihedralTypes detected: ", datafile.nDihedralTypes);
			fflush (stdout);
		}

		if (strstr (lineString, "improper types"))
		{
			sscanf (lineString, "%d \n", &datafile.nImproperTypes);
			fprintf(stdout, "%s%d\n", "nImproperTypes detected: ", datafile.nImproperTypes);
			fflush (stdout);
		}

		if ((datafile.nAtoms >= 0) && (datafile.nBonds >= 0) && (datafile.nAngles >= 0) && (datafile.nDihedrals >= 0) && (datafile.nImpropers >= 0) && (printHeaderInfo))
			printHeaderInfo = 0;

		if (strstr (lineString, "Atoms"))
		{
			isAtomLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Bonds"))
		{
			isBondLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Angles"))
		{
			isAngleLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Dihedrals"))
		{
			isDihedralLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (strstr (lineString, "Impropers"))
		{
			isImproperLine = 1;
			fgets (lineString, 1000, input);
			fgets (lineString, 1000, input);
		}

		if (isAtomLine)
		{
			sscanf (lineString, "%d %d %d %f %f %f %f\n", 
				&(*atoms)[nAtomLine].id, 
				&(*atoms)[nAtomLine].molType, 
				&(*atoms)[nAtomLine].atomType, 
				&(*atoms)[nAtomLine].charge, 
				&(*atoms)[nAtomLine].x, 
				&(*atoms)[nAtomLine].y, 
				&(*atoms)[nAtomLine].z);
			nAtomLine++;
			if (nAtomLine == datafile.nAtoms)
				isAtomLine = 0;
		}

		if (isBondLine)
		{
			sscanf (lineString, "%d %d %d %d\n", 
				&(*bonds)[nBondLine].id, 
				&(*bonds)[nBondLine].bondType, 
				&(*bonds)[nBondLine].atom1, 
				&(*bonds)[nBondLine].atom2);
			nBondLine++;
			if (nBondLine == datafile.nBonds)
				isBondLine = 0;
		}

		if (isAngleLine)
		{
			sscanf (lineString, "%d %d %d %d %d\n", 
				&(*angles)[nAngleLine].id, 
				&(*angles)[nAngleLine].angleType, 
				&(*angles)[nAngleLine].atom1, 
				&(*angles)[nAngleLine].atom2, 
				&(*angles)[nAngleLine].atom3);
			nAngleLine++;
			if (nAngleLine == datafile.nAngles)
				isAngleLine = 0;
		}

		if (isDihedralLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*dihedrals)[nDihedralLine].id, 
				&(*dihedrals)[nDihedralLine].dihedralType, 
				&(*dihedrals)[nDihedralLine].atom1, 
				&(*dihedrals)[nDihedralLine].atom2, 
				&(*dihedrals)[nDihedralLine].atom3, 
				&(*dihedrals)[nDihedralLine].atom4);
			nDihedralLine++;
			if (nDihedralLine == datafile.nDihedrals)
				isDihedralLine = 0;
		}

		if (isImproperLine)
		{
			sscanf (lineString, "%d %d %d %d %d %d\n", 
				&(*impropers)[nImproperLine].id, 
				&(*impropers)[nImproperLine].improperType, 
				&(*impropers)[nImproperLine].atom1, 
				&(*impropers)[nImproperLine].atom2, 
				&(*impropers)[nImproperLine].atom3, 
				&(*impropers)[nImproperLine].atom4);
			nImproperLine++;
			if (nImproperLine == datafile.nImpropers)
				isImproperLine = 0;
		}
	}

	printf("\nFrom input data file:\n\n nAtoms: %d\n nBonds: %d\n nAngles: %d\n nDihedrals: %d\n nImpropers: %d\n\n", datafile.nAtoms, datafile.nBonds, datafile.nAngles, datafile.nDihedrals, datafile.nImpropers);

	rewind (input);
	return datafile;
}

CONFIG *readConfig (FILE *inputConfigFile, int *nLines_return)
{
	rewind (inputConfigFile);
	CONFIG *inputVectors;
	char lineString[1000];
	int nLines = 0;

	while (fgets (lineString, 1000, inputConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			nLines++;
		}
	}

	inputVectors = (CONFIG *) malloc (nLines * sizeof (CONFIG));
	rewind (inputConfigFile);
	nLines = 0;

	while (fgets (lineString, 1000, inputConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			sscanf (lineString, "%d %d\n", &inputVectors[nLines].atom1, &inputVectors[nLines].atom2);
			nLines++;
		}
	}

	(*nLines_return) = nLines;

	rewind (inputConfigFile);
	return inputVectors;
}

CONFIG *readVWDRadius (FILE *inputVDWConfigFile, int *nLines_return)
{
	rewind (inputVDWConfigFile);
	CONFIG *inputVectors;
	char lineString[1000];
	int nLines = 0;

	while (fgets (lineString, 1000, inputVDWConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			nLines++;
		}
	}

	inputVectors = (CONFIG *) malloc (nLines * sizeof (CONFIG));
	rewind (inputVDWConfigFile);
	nLines = 0;

	while (fgets (lineString, 1000, inputVDWConfigFile) != NULL)
	{
		if (lineString[0] != '#')
		{
			sscanf (lineString, "%d %f\n", &inputVectors[nLines].atom1, &inputVectors[nLines].radius);
			nLines++;
		}
	}

	for (int i = 0; i < nLines; ++i)
	{
		inputVectors[i].radius /= 2;
	}

	(*nLines_return) = nLines;

	rewind (inputVDWConfigFile);
	return inputVectors;
}

DUMPFILE_INFO getDumpFileInfo (FILE *inputDumpFile)
{
	rewind (inputDumpFile);
	char lineString[1000];
	DUMPFILE_INFO dumpfile;

	for (int i = 0; i < 9; ++i)
	{
		fgets (lineString, 1000, inputDumpFile);
		if (i == 1)
		{
			sscanf (lineString, "%d", &dumpfile.timestep);
		}
		if (i == 3)
		{
			sscanf (lineString, "%d", &dumpfile.nAtoms);
		}
		if (i == 5)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.xlo, &dumpfile.xhi);
		}
		if (i == 6)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.ylo, &dumpfile.yhi);
		}
		if (i == 7)
		{
			sscanf (lineString, "%f %f\n", &dumpfile.zlo, &dumpfile.zhi);
		}
	}

	// Checking the information
	// printf("Timestep: %d\nNumber of atoms: %d\nxlo: %f; xhi: %f\nylo: %f; yhi: %f\nzlo: %f; zhi: %f\n", dumpfile.timestep, dumpfile.nAtoms, dumpfile.xlo, dumpfile.xhi, dumpfile.ylo, dumpfile.yhi, dumpfile.zlo, dumpfile.zhi);

	rewind (inputDumpFile);
	return dumpfile;
}

int countNTimeframes (FILE *inputDumpFile)
{
	printf("Counting the number of timeframes in input dump file...\n");
	fflush (stdout);

	rewind (inputDumpFile);

	int nLines = 0, nTimeframes = 0, nTimeframes2 = 0;
	char lineString[1000];

	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);
	int fault = 0, nAtoms, currentAtomID;

	while ((fgets (lineString, 1000, inputDumpFile) != NULL))
	{
		nLines++;

		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			startAgain:
			fault = 0;

			for (int i = 0; i < 8; ++i)
			{
				fgets (lineString, 1000, inputDumpFile);
				nLines++;

				if (i == 1 && strstr (lineString, "ITEM: NUMBER OF ATOMS") == 0)
					fault++;
				if (i == 2)
					sscanf (lineString, "%d\n", &nAtoms);
				if (i == 3 && strstr (lineString, "ITEM: BOX") == 0)
					fault++;
				if (i == 7 && strstr (lineString, "ITEM: ATOMS") == 0)
					fault++;
				if (strstr (lineString, "ITEM: TIMESTEP"))
					goto startAgain;
			}

			for (int i = 0; i < nAtoms; ++i)
			{
				fgets (lineString, 1000, inputDumpFile);
				nLines++;
				sscanf (lineString, "%d \n", &currentAtomID);

				if (currentAtomID != (i + 1))
					fault++;
				if (strstr (lineString, "ITEM: TIMESTEP"))
					goto startAgain;
			}

			// if fault is zero at the end of the timestep, then increment the variable
			if (fault == 0)
			{
				nTimeframes2++;
				if ((nTimeframes2 % 50) == 0)
				{
					printf("Scanning timestep: %5d    \r", nTimeframes2);
					fflush (stdout);
				}
			}
		}
	}

	rewind (inputDumpFile);

	printf("\nNumber of lines in the input dumpfile: %d\n", nLines);
	nTimeframes = nLines / (dumpfile.nAtoms + 9);
	printf("nTimeframes:  %d\n", nTimeframes);
	printf("nTimeframes2: %d\n", nTimeframes2);

	return nTimeframes2;
}