// Read dump file and data file
// Read the config file containing information about vectors
// Check the data file for bond connection.
// Draw a vector only if those two atoms are covalently connected.
// Save info about all the vectors
// Then calculate orientation order parameter

// This code should be written in a generic manner, such that it'll work for all systems

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

typedef struct nLines_config
{
	int nLines_inputVectors, nLines_vwdSize, nLines_freeVolumeconfig, nLines_HBondAtoms;
} NLINES_CONFIG;

typedef struct distVar
{
	float maxDist, binSize_OOP, binSize_deg, binSize_dist;
	float binStart_dist, binEnd_dist, binStart_OOP, binEnd_OOP, binStart_deg, binEnd_deg;
	int nBins_dist, nBins_OOP, nBins_deg, size_degrees, size_oop;

	// Optional element to add, for easy access
	int nElements;
} DIST_VAR;

typedef struct orderParameter
{
	int atom1, atom2, atom3, atom4;
	float distance, theta_rad, theta_deg, orderParameter;
} ORDERPARAMETER;

typedef struct dumpinfo
{
	int timestep, nAtoms;
	float xlo, xhi, ylo, yhi, zlo, zhi;
} DUMPFILE_INFO;

typedef struct config
{
	int atom1, atom2;
	float radius;
} CONFIG;

typedef struct datafile_atoms
{
	int resNumber, ix, iy, iz;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2, atom1Type, atom2Type;
	float x1, y1, z1, x2, y2, z2, xc, yc, zc;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

int isFile(const char *name)
{
	DIR *directory = opendir (name);
	if (directory!=NULL)
	{
		closedir(directory);
		return 0;
	}
	if(errno==ENOTDIR)
	{
		return 1;
	}

	return -1;
}

int displayFiles(const char *fileExtension)
{
	int nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;
	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			printf("%d --> %s\n", nFiles, filePointer -> d_name);
		}
	}
	return nFiles;
}

char *getInputFileName()
{
	int nFiles = 0;
	char *inputFileName, fileExtension[200];
	int fileRequired;
	inputFileName = (char *) malloc (200 * sizeof (char));

	getFilenameAgain:
	printf("Enter the file extension or a match string to search in current directory...\n --> "); scanf ("%s", fileExtension);
	printf("\n");
	nFiles = displayFiles (fileExtension);

	if (nFiles > 0)
	{
		fprintf(stdout, "\nWhich file would you like to input? Enter a number between (1 to %d): ", nFiles); 
		fflush (stdout);
		scanf ("%d", &fileRequired);
	}
	else
	{
		printf("No files found with the match string. Try again!\n"); goto getFilenameAgain;
	}

	nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;

	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			if (fileRequired == nFiles)
			{
				strcpy (inputFileName, filePointer -> d_name);
			}
		}
	}
	return inputFileName;
}

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

	rewind (inputDumpFile);
	return dumpfile;
}

float translatePeriodicDistance (float coord1, float coord2, float halfBoxDistance)
{
	if (abs (coord1 - coord2) > halfBoxDistance)
	{
		if (coord2 > coord1)
			coord2 -= (halfBoxDistance * 2);
		else if (coord1 > coord2)
			coord2 += (halfBoxDistance * 2);
	}
	return coord2;
}

float findBondCenter (float coord1, int image1, float coord2, int image2, float xlo, float xhi)
{
	float bondCenter;

	if (image1 > image2)
		coord1 = xhi + coord1 - xlo;
	else if (image1 < image2)
		coord1 = xlo - (xhi - coord1);

	bondCenter = ((coord1 + coord2) / 2);

	return bondCenter;
}

float findConnectedAtom_periodicTranslation (float coord1, int image1, float coord2, int image2, float xlo, float xhi)
{
	if (image1 > image2)
		coord1 = xhi + coord1 - xlo;
	else if (image1 < image2)
		coord1 = xlo - (xhi - coord1);

	return coord1;
}

ORDERPARAMETER *computeOrderParameter (DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, int currentTimestep, unsigned int nElements)
{
	FILE *allData;
	char *allData_string;
	allData_string = (char *) malloc (50 * sizeof (char));
	sprintf (allData_string, "logs/allData_%d.oop", currentTimestep);
	allData = fopen (allData_string, "w");

	fprintf(allData, "atom1, atom2, atom3, atom4, distance, angle (rad), angle (deg), OOP\n");

	ORDERPARAMETER *allData_array;
	allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));

	unsigned int currentElement = 0;

	float x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, distance, dotProduct, magnitude1, magnitude2, cosTheta, theta, orderParameter;
	float xDistHalf = ((dumpfile.xhi - dumpfile.xlo) / 2), yDistHalf = ((dumpfile.yhi - dumpfile.ylo) / 2), zDistHalf = ((dumpfile.zhi - dumpfile.zlo) / 2);

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		bonds[i].atom1Type = dumpAtoms[bonds[i].atom1 - 1].atomType;
		bonds[i].atom2Type = dumpAtoms[bonds[i].atom2 - 1].atomType;

		// Checking if the bonds correspond to inputVectors[0]; from the first line of the config file
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom1Type == inputVectors[0].atom2 && bonds[i].atom2Type == inputVectors[0].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom1Type == inputVectors[1].atom2 && bonds[j].atom2Type == inputVectors[1].atom1))
				{
					// Finding the center of two bonds (x1, y1, z1) and (x2, y2, z2)
					// x1 = (dumpAtoms[bonds[i].atom1 - 1].x + dumpAtoms[bonds[i].atom2 - 1].x) / 2; 
					// y1 = (dumpAtoms[bonds[i].atom1 - 1].y + dumpAtoms[bonds[i].atom2 - 1].y) / 2; 
					// z1 = (dumpAtoms[bonds[i].atom1 - 1].z + dumpAtoms[bonds[i].atom2 - 1].z) / 2; 

					x1 = findBondCenter (dumpAtoms[bonds[i].atom1 - 1].x, dumpAtoms[bonds[i].atom1 - 1].ix, dumpAtoms[bonds[i].atom2 - 1].x, dumpAtoms[bonds[i].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
					y1 = findBondCenter (dumpAtoms[bonds[i].atom1 - 1].y, dumpAtoms[bonds[i].atom1 - 1].iy, dumpAtoms[bonds[i].atom2 - 1].y, dumpAtoms[bonds[i].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
					z1 = findBondCenter (dumpAtoms[bonds[i].atom1 - 1].z, dumpAtoms[bonds[i].atom1 - 1].iz, dumpAtoms[bonds[i].atom2 - 1].z, dumpAtoms[bonds[i].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

					// x2 = (dumpAtoms[bonds[j].atom1 - 1].x + dumpAtoms[bonds[j].atom2 - 1].x) / 2; 
					// y2 = (dumpAtoms[bonds[j].atom1 - 1].y + dumpAtoms[bonds[j].atom2 - 1].y) / 2; 
					// z2 = (dumpAtoms[bonds[j].atom1 - 1].z + dumpAtoms[bonds[j].atom2 - 1].z) / 2;

					x2 = findBondCenter (dumpAtoms[bonds[j].atom1 - 1].x, dumpAtoms[bonds[j].atom1 - 1].ix, dumpAtoms[bonds[j].atom2 - 1].x, dumpAtoms[bonds[j].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
					y2 = findBondCenter (dumpAtoms[bonds[j].atom1 - 1].y, dumpAtoms[bonds[j].atom1 - 1].iy, dumpAtoms[bonds[j].atom2 - 1].y, dumpAtoms[bonds[j].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
					z2 = findBondCenter (dumpAtoms[bonds[j].atom1 - 1].z, dumpAtoms[bonds[j].atom1 - 1].iz, dumpAtoms[bonds[j].atom2 - 1].z, dumpAtoms[bonds[j].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

					x2 = translatePeriodicDistance (x1, x2, xDistHalf); 
					y2 = translatePeriodicDistance (y1, y2, yDistHalf); 
					z2 = translatePeriodicDistance (z1, z2, zDistHalf);

					// Distance between the centers of two bonds
					distance = sqrt (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2  - z1)));

					// Storing the positions of all 4 atoms forming the two bonds of interest
					x1 = dumpAtoms[bonds[i].atom1 - 1].x; 
					y1 = dumpAtoms[bonds[i].atom1 - 1].y; 
					z1 = dumpAtoms[bonds[i].atom1 - 1].z;

					x2 = dumpAtoms[bonds[i].atom2 - 1].x; 
					y2 = dumpAtoms[bonds[i].atom2 - 1].y; 
					z2 = dumpAtoms[bonds[i].atom2 - 1].z; 

					x3 = dumpAtoms[bonds[j].atom1 - 1].x; 
					y3 = dumpAtoms[bonds[j].atom1 - 1].y; 
					z3 = dumpAtoms[bonds[j].atom1 - 1].z; 

					x4 = dumpAtoms[bonds[j].atom2 - 1].x; 
					y4 = dumpAtoms[bonds[j].atom2 - 1].y; 
					z4 = dumpAtoms[bonds[j].atom2 - 1].z; 

					dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3)) + ((z2 - z1) * (z4 - z3)); 
					magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)); 
					magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)) + ((z4 - z3) * (z4 - z3)); 

					cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2)); 
					theta = acosf (cosTheta); 
					orderParameter = ((3.0 * cosTheta * cosTheta) - 1.0) / 2.0;

					fprintf(allData, "%d %d %d %d %f %f %f %f\n", bonds[i].atom1, bonds[i].atom2, bonds[j].atom1, bonds[j].atom2, distance, theta, theta * 57.2958, orderParameter);

					allData_array[currentElement].atom1 = bonds[i].atom1; 
					allData_array[currentElement].atom2 = bonds[i].atom2; 
					allData_array[currentElement].atom3 = bonds[j].atom1; 
					allData_array[currentElement].atom4 = bonds[j].atom2; 
					allData_array[currentElement].distance = distance; 
					allData_array[currentElement].theta_rad = theta; 
					allData_array[currentElement].theta_deg = theta * 57.2958; 
					allData_array[currentElement].orderParameter = orderParameter; 

					currentElement++;
				}
			}
		}

		// Checking if the bonds correspond to inputVectors[1]; from the second line of the config file
		else if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom1Type == inputVectors[1].atom2 && bonds[i].atom2Type == inputVectors[1].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[0].atom1 && bonds[j].atom2Type == inputVectors[0].atom2) || (bonds[j].atom1Type == inputVectors[0].atom2 && bonds[j].atom2Type == inputVectors[0].atom1))
				{
					// Finding the center of two bonds (x1, y1, z1) and (x2, y2, z2)
					x1 = (dumpAtoms[bonds[i].atom1 - 1].x + dumpAtoms[bonds[i].atom2 - 1].x) / 2; y1 = (dumpAtoms[bonds[i].atom1 - 1].y + dumpAtoms[bonds[i].atom2 - 1].y) / 2; z1 = (dumpAtoms[bonds[i].atom1 - 1].z + dumpAtoms[bonds[i].atom2 - 1].z) / 2; 

					x2 = (dumpAtoms[bonds[j].atom1 - 1].x + dumpAtoms[bonds[j].atom2 - 1].x) / 2; y2 = (dumpAtoms[bonds[j].atom1 - 1].y + dumpAtoms[bonds[j].atom2 - 1].y) / 2; z2 = (dumpAtoms[bonds[j].atom1 - 1].z + dumpAtoms[bonds[j].atom2 - 1].z) / 2;

					x2 = translatePeriodicDistance (x1, x2, xDistHalf); y2 = translatePeriodicDistance (y1, y2, yDistHalf); z2 = translatePeriodicDistance (z1, z2, zDistHalf);

					// Distance between the centers of two bonds
					distance = sqrt (((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2  - z1)));

					// Storing the positions of all 4 atoms forming the two bonds of interest
					x1 = dumpAtoms[bonds[i].atom1 - 1].x; 
					y1 = dumpAtoms[bonds[i].atom1 - 1].y; 
					z1 = dumpAtoms[bonds[i].atom1 - 1].z; 
					x2 = dumpAtoms[bonds[i].atom2 - 1].x; 
					y2 = dumpAtoms[bonds[i].atom2 - 1].y; 
					z2 = dumpAtoms[bonds[i].atom2 - 1].z; 
					x3 = dumpAtoms[bonds[j].atom1 - 1].x; 
					y3 = dumpAtoms[bonds[j].atom1 - 1].y; 
					z3 = dumpAtoms[bonds[j].atom1 - 1].z; 
					x4 = dumpAtoms[bonds[j].atom2 - 1].x; 
					y4 = dumpAtoms[bonds[j].atom2 - 1].y; 
					z4 = dumpAtoms[bonds[j].atom2 - 1].z; 

					dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3)) + ((z2 - z1) * (z4 - z3)); 
					magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)); 
					magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)) + ((z4 - z3) * (z4 - z3)); 
					cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2)); theta = acosf (cosTheta); 
					orderParameter = ((3 * cosTheta * cosTheta) - 1) / 2;

					// fprintf(allData, "%d %d %d %d %f %f %f %f\n", bonds[i].atom1, bonds[i].atom2, bonds[j].atom1, bonds[j].atom2, distance, theta, theta * 57.2958, orderParameter);

					allData_array[currentElement].atom1 = bonds[i].atom1; 
					allData_array[currentElement].atom2 = bonds[i].atom2; 
					allData_array[currentElement].atom3 = bonds[j].atom1; 
					allData_array[currentElement].atom4 = bonds[j].atom2; 
					allData_array[currentElement].distance = distance; 
					allData_array[currentElement].theta_rad = theta; 
					allData_array[currentElement].theta_deg = theta * 57.2958; 
					allData_array[currentElement].orderParameter = orderParameter; 
					currentElement++;
				}				
			}
		}
	}

	// fclose (allData);
	return allData_array;
}

unsigned int getNElements (DATAFILE_INFO datafile, DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, CONFIG *inputVectors)
{
	unsigned int nElements = 0;

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType;
		bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;

		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom1Type == inputVectors[0].atom2 && bonds[i].atom2Type == inputVectors[0].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom1Type == inputVectors[1].atom2 && bonds[j].atom2Type == inputVectors[1].atom1))
				{
					nElements++;
				}
			}
		}
		else if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom1Type == inputVectors[1].atom2 && bonds[i].atom2Type == inputVectors[1].atom1))
		{
			for (int j = (i + 1); j < datafile.nBonds; ++j)
			{
				bonds[j].atom1Type = dumpAtoms[bonds[j].atom1 - 1].atomType;
				bonds[j].atom2Type = dumpAtoms[bonds[j].atom2 - 1].atomType;

				if ((bonds[j].atom1Type == inputVectors[0].atom1 && bonds[j].atom2Type == inputVectors[0].atom2) || (bonds[j].atom1Type == inputVectors[0].atom2 && bonds[j].atom2Type == inputVectors[0].atom1))
				{
					nElements++;
				}				
			}
		}

	}

	return nElements;
}

typedef struct distribution
{
	float binStart_OOP, binEnd_OOP, binStart_dist, binEnd_dist, binStart_deg, binEnd_deg;
	int count;
} DISTRIBUTION;

void setDistributionZero (DISTRIBUTION **rawArray, int arraySize)
{
	for (int i = 0; i < arraySize; ++i)
	{
		(*rawArray)[i].count = 0; (*rawArray)[i].binStart_OOP = 0; (*rawArray)[i].binEnd_OOP = 0; (*rawArray)[i].binStart_dist = 0; (*rawArray)[i].binEnd_dist = 0; (*rawArray)[i].binStart_deg = 0; (*rawArray)[i].binEnd_deg = 0;
	}
}

int getIndex1d (int i, int j, int width)
{
	// width is the max value of 'j'
	int index1d = 0;
	index1d = (width * i) + j;
	return index1d;
}

long long int getIndex1d_from3d (int x, int xWidth, int y, int yWidth, int z, int zWidth)
{
	long long int arrayIndex;
	arrayIndex = (long long int) ((z * yWidth * xWidth) + (y * xWidth) + x);
	return arrayIndex;
}

void computeDistribution_OOP (ORDERPARAMETER *allData_array, DIST_VAR plotVars, DISTRIBUTION **distribution_OOP, int nThreads)
{
	DIST_VAR currentBounds;
	currentBounds.binStart_OOP = plotVars.binStart_OOP;
	currentBounds.binStart_dist = plotVars.binStart_dist;

	// OMP setup
	omp_set_num_threads(nThreads);

	int index1d;

	for (int i = 0; i < plotVars.nBins_OOP; ++i)
	{
		fprintf(stdout, "Computing OOP distribution... %d/%d            \r", i, plotVars.nBins_OOP);
		fflush (stdout);

		currentBounds.binEnd_OOP = currentBounds.binStart_OOP + plotVars.binSize_OOP;
		currentBounds.binStart_dist = plotVars.binStart_dist;

		for (int j = 0; j < plotVars.nBins_dist; ++j)
		{
			currentBounds.binEnd_dist = currentBounds.binStart_dist + plotVars.binSize_dist;

			#pragma omp parallel for
			for (int k = 0; k < plotVars.nElements; ++k)
			{
				if (allData_array[k].orderParameter <= currentBounds.binEnd_OOP && allData_array[k].orderParameter > currentBounds.binStart_OOP && allData_array[k].distance <= currentBounds.binEnd_dist && allData_array[k].distance > currentBounds.binStart_dist)
				{
					index1d = getIndex1d (i, j, plotVars.nBins_dist);

					(*distribution_OOP)[index1d].count++;
					(*distribution_OOP)[index1d].binStart_OOP = currentBounds.binStart_OOP;
					(*distribution_OOP)[index1d].binEnd_OOP = currentBounds.binEnd_OOP;
					(*distribution_OOP)[index1d].binStart_dist = currentBounds.binStart_dist;
					(*distribution_OOP)[index1d].binEnd_dist = currentBounds.binEnd_dist;
					(*distribution_OOP)[index1d].binStart_deg = 0;
					(*distribution_OOP)[index1d].binEnd_deg = 0;
				}
			}
			currentBounds.binStart_dist = currentBounds.binEnd_dist;
		}
		currentBounds.binStart_OOP = currentBounds.binEnd_OOP;
	}
}

void computeDistribution_theta (ORDERPARAMETER *allData_array, DIST_VAR plotVars, DISTRIBUTION **distribution_degrees, int nThreads)
{
	DIST_VAR currentBounds;
	currentBounds.binStart_deg = plotVars.binStart_deg;
	currentBounds.binStart_dist = plotVars.binStart_dist;

	// OMP setup
	omp_set_num_threads(nThreads);

	int index1d;

	for (int i = 0; i < plotVars.nBins_deg; ++i)
	{
		fprintf(stdout, "Computing theta distribution... %d/%d               \r", i, plotVars.nBins_deg);
		fflush (stdout);

		currentBounds.binEnd_deg = currentBounds.binStart_deg + plotVars.binSize_deg;
		currentBounds.binStart_dist = plotVars.binStart_dist;

		for (int j = 0; j < plotVars.nBins_dist; ++j)
		{
			currentBounds.binEnd_dist = currentBounds.binStart_dist + plotVars.binSize_dist;

			#pragma omp parallel for
			for (int k = 0; k < plotVars.nElements; ++k)
			{
				if (allData_array[k].theta_deg <= currentBounds.binEnd_deg && allData_array[k].theta_deg > currentBounds.binStart_deg && allData_array[k].distance <= currentBounds.binEnd_dist && allData_array[k].distance > currentBounds.binStart_dist)
				{
					index1d = getIndex1d (i, j, plotVars.nBins_dist);

					(*distribution_degrees)[index1d].count++;
					(*distribution_degrees)[index1d].binStart_deg = currentBounds.binStart_deg;
					(*distribution_degrees)[index1d].binEnd_deg = currentBounds.binEnd_deg;
					(*distribution_degrees)[index1d].binStart_dist = currentBounds.binStart_dist;
					(*distribution_degrees)[index1d].binEnd_dist = currentBounds.binEnd_dist;
					(*distribution_degrees)[index1d].binStart_deg = 0;
					(*distribution_degrees)[index1d].binEnd_deg = 0;
				}
			}
			currentBounds.binStart_dist = currentBounds.binEnd_dist;
		}
		currentBounds.binStart_deg = currentBounds.binEnd_deg;
	}
}

void printDistribution_OOP (DISTRIBUTION *distribution_OOP, DIST_VAR plotVars)
{
	FILE *file_distribution_OOP, *file_distribution_OOP_info;
	file_distribution_OOP = fopen ("orderParameter.dist", "w");
	file_distribution_OOP_info = fopen ("orderParameter.dist.info", "w");

	int index1d, oop_index, dist_index;

	// Printing header information
	fprintf(file_distribution_OOP_info, "binStart_dist: %f\nbinEnd_dist: %f\nbinStart_OOP: %f\nbinEnd_OOP: %f\nnBins_dist: %d\nnBins_OOP: %d\nsize_oop: %d\n\n",
		plotVars.binStart_dist,
		plotVars.binEnd_dist,
		plotVars.binStart_OOP,
		plotVars.binEnd_OOP,
		plotVars.nBins_dist,
		plotVars.nBins_OOP,
		plotVars.size_oop);

	for (int oop_index = 0; oop_index < plotVars.nBins_OOP; ++oop_index)
	{
		fprintf(file_distribution_OOP, "\n");
		for (int dist_index = 0; dist_index < plotVars.nBins_dist; ++dist_index)
		{
			index1d = getIndex1d (oop_index, dist_index, plotVars.nBins_dist);
			fprintf(file_distribution_OOP, "%d\t", distribution_OOP[index1d].count);
		}
	}

	fclose (file_distribution_OOP);
	fclose (file_distribution_OOP_info);
}

void printDistribution_degrees (DISTRIBUTION *distribution_degrees, DIST_VAR plotVars)
{
	FILE *file_distribution_degrees, *file_distribution_degrees_info;
	file_distribution_degrees = fopen ("degrees.dist", "w");
	file_distribution_degrees_info = fopen ("degrees.dist.info", "w");

	// Printing header information
	fprintf(file_distribution_degrees_info, "binStart_dist: %f\nbinEnd_dist: %f\nbinStart_deg: %f\nbinEnd_deg: %f\nnBins_dist: %d\nnBins_deg: %d\nsize_degrees: %f\n\n",
		plotVars.binStart_dist,
		plotVars.binEnd_dist,
		plotVars.binStart_deg,
		plotVars.binEnd_deg,
		plotVars.nBins_dist,
		plotVars.nBins_deg,
		plotVars.size_degrees);
	
	int index1d, deg_index, dist_index;

	for (int deg_index = 0; deg_index < plotVars.nBins_deg; ++deg_index)
	{
		fprintf(file_distribution_degrees, "\n");
		for (int dist_index = 0; dist_index < plotVars.nBins_dist; ++dist_index)
		{
			index1d = getIndex1d (deg_index, dist_index, plotVars.nBins_dist);
			fprintf(file_distribution_degrees, "%d\t", distribution_degrees[index1d].count);
		}
	}

	fclose (file_distribution_degrees);
	fclose (file_distribution_degrees_info);

	// Printing normalized data
	FILE *file_distribution_degrees_normalized;
	file_distribution_degrees_normalized = fopen ("degrees.dist.norm", "w");

	int *maxCount;
	maxCount = (int *) calloc (plotVars.nBins_dist, sizeof (int));

	// Finding the max count in every distance bin
	for (int dist_index = 0; dist_index < plotVars.nBins_dist; ++dist_index)
	{
		for (int deg_index = 0; deg_index < plotVars.nBins_deg; ++deg_index)
		{
			index1d = getIndex1d (deg_index, dist_index, plotVars.nBins_dist);
			if (distribution_degrees[index1d].count > maxCount[dist_index])
				maxCount[dist_index] = distribution_degrees[index1d].count;
		}
	}

	// Printing the normalized output values
	float normalizedDistribution;

	for (int deg_index = 0; deg_index < plotVars.nBins_deg; ++deg_index)
	{
		fprintf(file_distribution_degrees_normalized, "\n");
		for (int dist_index = 0; dist_index < plotVars.nBins_dist; ++dist_index)
		{
			index1d = getIndex1d (deg_index, dist_index, plotVars.nBins_dist);
			if (maxCount[dist_index] > 0)
				normalizedDistribution = (float) distribution_degrees[index1d].count / (float) maxCount[dist_index];
			else
				normalizedDistribution = 0;

			fprintf(file_distribution_degrees_normalized, "%f\t", normalizedDistribution);
		}
	}

	fclose (file_distribution_degrees_normalized);
}

void computeBondRDF (DATA_ATOMS *dumpAtoms, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, DATA_BONDS *bonds, CONFIG *inputVectors, DIST_VAR plotVars, int nThreads, float binSize_dist_RDF, float **bondRDF, int *RDFcounter, int currentTimestep)
{
	FILE *bondRDF_logfile;
	char *bondRDF_logfilename;
	bondRDF_logfilename = (char *) malloc (50 * sizeof (char));
	snprintf (bondRDF_logfilename, 50, "bondRDF_logs/%d.rdf", currentTimestep);
	bondRDF_logfile = fopen (bondRDF_logfilename, "w");

	(*RDFcounter)++;

	float xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo), simVolume = (xDist * yDist * zDist), bondDensity;
	float xDistHalf = (xDist / 2), yDistHalf = (yDist / 2), zDistHalf = (zDist / 2);
	int nBonds_RDF = 0;

	// Computing bondRDF
	float binStart_dist_RDF = 0, binEnd_dist_RDF, distance;
	int nBins_dist_RDF = (int) (plotVars.maxDist / binSize_dist_RDF);
	int *nBonds_inBin_dist_RDF;
	nBonds_inBin_dist_RDF = (int *) calloc (nBins_dist_RDF, sizeof (int));
	float *nBonds_inBin_dist_RDF_float;
	nBonds_inBin_dist_RDF_float = (float *) malloc (nBins_dist_RDF * sizeof (float));

	float x_translated, y_translated, z_translated;

	omp_set_num_threads (nThreads);
	#pragma omp parallel for
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom2Type == inputVectors[0].atom1 && bonds[i].atom1Type == inputVectors[0].atom2))
		{
			bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType; bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;
			bonds[i].x1 = dumpAtoms[bonds[i].atom1 - 1].x; bonds[i].y1 = dumpAtoms[bonds[i].atom1 - 1].y; bonds[i].z1 = dumpAtoms[bonds[i].atom1 - 1].z;
			bonds[i].x2 = dumpAtoms[bonds[i].atom2 - 1].x; bonds[i].y2 = dumpAtoms[bonds[i].atom2 - 1].y; bonds[i].z2 = dumpAtoms[bonds[i].atom2 - 1].z;

			bonds[i].xc = findBondCenter (bonds[i].x1, dumpAtoms[bonds[i].atom1 - 1].ix, bonds[i].x2, dumpAtoms[bonds[i].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
			bonds[i].yc = findBondCenter (bonds[i].y1, dumpAtoms[bonds[i].atom1 - 1].iy, bonds[i].y2, dumpAtoms[bonds[i].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
			bonds[i].zc = findBondCenter (bonds[i].z1, dumpAtoms[bonds[i].atom1 - 1].iz, bonds[i].z2, dumpAtoms[bonds[i].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

			// bonds[i].xc = (bonds[i].x1 + bonds[i].x2) / 2; bonds[i].yc = (bonds[i].y1 + bonds[i].y2) / 2; bonds[i].zc = (bonds[i].z1 + bonds[i].z2) / 2;

			for (int j = 0; j < datafile.nBonds; ++j)
			{
				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom2Type == inputVectors[1].atom1 && bonds[j].atom1Type == inputVectors[1].atom2))
				{
					bonds[j].atom1Type = (int) dumpAtoms[bonds[j].atom1 - 1].atomType; bonds[j].atom2Type = (int) dumpAtoms[bonds[j].atom2 - 1].atomType;
					bonds[j].x1 = dumpAtoms[bonds[j].atom1 - 1].x; bonds[j].y1 = dumpAtoms[bonds[j].atom1 - 1].y; bonds[j].z1 = dumpAtoms[bonds[j].atom1 - 1].z;
					bonds[j].x2 = dumpAtoms[bonds[j].atom2 - 1].x; bonds[j].y2 = dumpAtoms[bonds[j].atom2 - 1].y; bonds[j].z2 = dumpAtoms[bonds[j].atom2 - 1].z;

					bonds[j].xc = findBondCenter (bonds[j].x1, dumpAtoms[bonds[j].atom1 - 1].ix, bonds[j].x2, dumpAtoms[bonds[j].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
					bonds[j].yc = findBondCenter (bonds[j].y1, dumpAtoms[bonds[j].atom1 - 1].iy, bonds[j].y2, dumpAtoms[bonds[j].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
					bonds[j].zc = findBondCenter (bonds[j].z1, dumpAtoms[bonds[j].atom1 - 1].iz, bonds[j].z2, dumpAtoms[bonds[j].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

					// bonds[j].xc = (bonds[j].x1 + bonds[j].x2) / 2; bonds[j].yc = (bonds[j].y1 + bonds[j].y2) / 2; bonds[j].zc = (bonds[j].z1 + bonds[j].z2) / 2;

					x_translated = translatePeriodicDistance (bonds[i].xc, bonds[j].xc, xDistHalf);
					y_translated = translatePeriodicDistance (bonds[i].yc, bonds[j].yc, yDistHalf);
					z_translated = translatePeriodicDistance (bonds[i].zc, bonds[j].zc, zDistHalf);

					// if ((x_translated != bonds[j].xc) || (y_translated != bonds[j].yc) || (z_translated != bonds[j].zc))
					// {
					// 	printf("before translation => [%.3f, %.3f, %.3f] <=> ref: [%.3f, %.3f, %.3f]\nafter translation => [%.3f, %.3f, %.3f] <=> ref: [%.3f, %.3f, %.3f]\n\n", bonds[j].xc, bonds[j].yc, bonds[j].zc, bonds[i].xc, bonds[i].yc, bonds[i].zc, x_translated, y_translated, z_translated, bonds[i].xc, bonds[i].yc, bonds[i].zc);
					// 	sleep (1);
					// }

					// distance = sqrt (pow ((bonds[i].xc - bonds[j].xc), 2) + pow ((bonds[i].yc - bonds[j].yc), 2) + pow ((bonds[i].zc - bonds[j].zc), 2));
					distance = sqrt (pow ((bonds[i].xc - x_translated), 2) + pow ((bonds[i].yc - y_translated), 2) + pow ((bonds[i].zc - z_translated), 2));
					binStart_dist_RDF = 0.0;

					fprintf(bondRDF_logfile, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", bonds[i].xc, bonds[i].yc, bonds[i].zc, bonds[j].xc, bonds[j].yc, bonds[j].zc, x_translated, y_translated, z_translated, distance);

					for (int k = 0; k < nBins_dist_RDF; ++k)
					{
						binEnd_dist_RDF = binStart_dist_RDF + binSize_dist_RDF;

						// If i = j, then distance = 0.0
						// neglect the bond pairs where the distance = 0.0
						if (distance > 0 && distance <= binEnd_dist_RDF && distance > binStart_dist_RDF)
							nBonds_inBin_dist_RDF[k]++;

						binStart_dist_RDF = binEnd_dist_RDF;
					}
				}
			}
		}
	}

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		// Finding the bulk density of bond present in inputVectors[1]
		if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom2Type == inputVectors[1].atom1 && bonds[i].atom1Type == inputVectors[1].atom2))
			nBonds_RDF++;		
	}

	bondDensity = (float) nBonds_RDF * 20.0 / simVolume;

	binStart_dist_RDF = 0.0;

	for (int i = 0; i < nBins_dist_RDF; ++i)
	{
		binEnd_dist_RDF = binStart_dist_RDF + binSize_dist_RDF;
		nBonds_inBin_dist_RDF_float[i] = (float) nBonds_inBin_dist_RDF[i] / (4.0 * 3.14 * pow (binStart_dist_RDF, 2) * (binEnd_dist_RDF - binStart_dist_RDF));
		(*bondRDF)[i] += (float) nBonds_inBin_dist_RDF_float[i] / bondDensity;
		binStart_dist_RDF = binEnd_dist_RDF;
	}

	fclose (bondRDF_logfile);
}

void printBondRDF (float *bondRDF, int RDFcounter, int nBins_dist_RDF, float binSize_dist_RDF)
{
	FILE *file_bondRDF;
	file_bondRDF = fopen ("bondRDF.output", "w");

	for (int i = 0; i < nBins_dist_RDF; ++i)
	{
		fprintf(file_bondRDF, "%f %f\n", 
			binSize_dist_RDF * (float) i,
			bondRDF[i]/RDFcounter);
	}

	fclose (file_bondRDF);
}

typedef struct freeVolumeVars
{
	float minProbeSize, maxProbeSize, delProbeSize, currentProbeSize;
	int nBins_probeSweep;
	float xLength, yLength, zLength;
	int nBins_dist_x, nBins_dist_y, nBins_dist_z;
	float delDistance;

	float binStart_dist, binEnd_dist, binSize_dist;
	int nBins_dist;
} FREEVOLUME_VARS;

FREEVOLUME_VARS getFreeVolumeVars (DUMPFILE_INFO dumpfile)
{
	FREEVOLUME_VARS freeVolumeVars;

	// Variables regarding probe size
	printf("Enter the minimum size of probe:\t");
	scanf ("%f", &freeVolumeVars.minProbeSize);
	printf("Enter the maximum size of probe:\t");
	scanf ("%f", &freeVolumeVars.maxProbeSize);
	printf("Enter the del (size) for probe:\t");
	scanf ("%f", &freeVolumeVars.delProbeSize);

	freeVolumeVars.nBins_probeSweep = (int) ((freeVolumeVars.maxProbeSize - freeVolumeVars.minProbeSize) / freeVolumeVars.delProbeSize) + 1;

	// Variables regarding 3D probe movement
	// delDistance can usually be set at 1/4th of probe radius
	freeVolumeVars.xLength = (dumpfile.xhi - dumpfile.xlo); freeVolumeVars.yLength = (dumpfile.yhi - dumpfile.ylo); freeVolumeVars.zLength = (dumpfile.zhi - dumpfile.zlo);

	// Variables for calculating free volume distribution
	freeVolumeVars.binStart_dist = 0;
	printf("Free volume distribution will be calculated from the center of respective atoms\n\n");
	printf("Enter the max. distance to consider for free volume distribution (float type):\t"); scanf ("%f", &freeVolumeVars.binEnd_dist);
	printf("Enter the bin size (number of bins will be calculated based on max. distance and bin size):\t"); scanf ("%f", &freeVolumeVars.binSize_dist);
	freeVolumeVars.nBins_dist = (int) ((freeVolumeVars.binEnd_dist - freeVolumeVars.binStart_dist) / freeVolumeVars.binSize_dist) + 1;

	return freeVolumeVars;
}

typedef struct freeVolumeDistribution
{
	float binStart_dist, binEnd_dist;
	int nUnoccupied, nOccupied, atomType;
} FREEVOLUME_DISTRIBUTION;

// i corresponds to probeSize
int *computeFreeVolume_checkOccupation (int i, DUMPFILE_INFO dumpfile, DATA_ATOMS *dumpAtoms, FREEVOLUME_VARS freeVolumeVars, CONFIG *vwdSize, int nThreads)
{
	DATA_ATOMS probePosition;
	float distance, x1, y1, z1, x2, y2, z2;
	int *isOccupied;
	long long int arraySize = (long long int)(freeVolumeVars.nBins_dist_x + 1) * (long long int)(freeVolumeVars.nBins_dist_y + 1) * (long long int)(freeVolumeVars.nBins_dist_z + 1), index1d, progress = 0;
	float progressPercent;
	isOccupied = (int *) calloc (arraySize, sizeof (int));

	omp_set_num_threads (nThreads);

	probePosition.x = dumpfile.xlo;

	float xDistHalf = (dumpfile.xhi - dumpfile.xlo) / 2, yDistHalf = (dumpfile.yhi - dumpfile.ylo) / 2, zDistHalf = (dumpfile.zhi - dumpfile.zlo) / 2;

	#pragma omp parallel for
	for (int j = 0; j < freeVolumeVars.nBins_dist_x; ++j)
	{
		probePosition.y = dumpfile.ylo;
		for (int k = 0; k < freeVolumeVars.nBins_dist_y; ++k)
		{
			probePosition.z = dumpfile.zlo;
			for (int l = 0; l < freeVolumeVars.nBins_dist_z; ++l)
			{
				progress++;
				progressPercent = (float)progress / (float)arraySize;
				printf("Checking site occupation for probeSize: %.2f... %3.4f %% completed            \r\r", freeVolumeVars.currentProbeSize, progressPercent * 100.0);
				fflush (stdout);

				for (int m = 0; m < dumpfile.nAtoms; ++m)
				{
					x1 = dumpAtoms[m].x; y1 = dumpAtoms[m].y; z1 = dumpAtoms[m].z;
					x2 = probePosition.x; y2 = probePosition.y; z2 = probePosition.z;

					x2 = translatePeriodicDistance (x1, x2, xDistHalf);
					y2 = translatePeriodicDistance (y1, y2, yDistHalf);
					z2 = translatePeriodicDistance (z1, z2, zDistHalf);
					
					distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

					index1d = getIndex1d_from3d (j, freeVolumeVars.nBins_dist_x, k, freeVolumeVars.nBins_dist_y, l, freeVolumeVars.nBins_dist_z);
					// printf("Checking occupation => index1d: %lld/%lld;", index1d, arraySize);
					// fflush (stdout);
					// printf(" r: %.2f\n", vwdSize[dumpAtoms[m].atomType - 1].radius);
					// fflush (stdout);

					if (index1d >= arraySize)
					{
						printf("Checking occupation => index1d: %lld/%lld;", index1d, arraySize);
						printf(" r: %.2f\n", vwdSize[dumpAtoms[m].atomType - 1].radius);
						arraySize += 100;
						isOccupied = (int *) realloc (isOccupied, sizeof (int) * arraySize);
					}

					if (distance < (freeVolumeVars.currentProbeSize + vwdSize[dumpAtoms[m].atomType - 1].radius))
					{
						isOccupied[index1d] = 1;
					}
				}
				probePosition.z += freeVolumeVars.delDistance;
			}
			probePosition.y += freeVolumeVars.delDistance;
		}
		probePosition.x += freeVolumeVars.delDistance;
	}

	printf("\n");

	return isOccupied;
}

// i corresponds to probeSize and j corresponds to atom type provided in the input config file
void computeFreeVolume_getDistribution (int i, int j, FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars, int *isOccupied, DUMPFILE_INFO dumpfile, DATA_ATOMS *dumpAtoms, CONFIG *freeVolumeconfig, int nThreads)
{
	DATA_ATOMS probePosition;
	float distance, x1, y1, z1, x2, y2, z2;
	// int *isOccupied;
	long long int index1d, arraySize = (long long int)freeVolumeVars.nBins_dist_x * (long long int)freeVolumeVars.nBins_dist_y * (long long int)freeVolumeVars.nBins_dist_z, progress = 0;
	float progressPercent;
	// isOccupied = (int *) malloc (arraySize * sizeof (int));

	omp_set_num_threads (nThreads);

	#pragma omp parallel for
	for (int k = 0; k < dumpfile.nAtoms; ++k)
	{
		progress++;
		printf("Computing distribution (probeSize: %.2f): %3.4f %% completed...           \r", freeVolumeVars.currentProbeSize, ((float) progress / (float) dumpfile.nAtoms) * 100.0);
		fflush (stdout);
		// If the atom type matches the atom type given in config file, then proceed with the calculation
		if (dumpAtoms[k].atomType == freeVolumeconfig[j].atom1)
		{
			x1 = dumpAtoms[k].x; y1 = dumpAtoms[k].y; z1 = dumpAtoms[k].z;

			// Going through all positions of the probe
			probePosition.x = dumpfile.xlo;
			for (int l = 0; l < freeVolumeVars.nBins_dist_x; ++l)
			{
				probePosition.y = dumpfile.ylo;
				for (int m = 0; m < freeVolumeVars.nBins_dist_y; ++m)
				{
					probePosition.z = dumpfile.zlo;
					for (int n = 0; n < freeVolumeVars.nBins_dist_z; ++n)
					{
						// Check if the position is unoccupied. Proceed only if it is unoccupied
						index1d = getIndex1d_from3d (l, freeVolumeVars.nBins_dist_x, m, freeVolumeVars.nBins_dist_y, n, freeVolumeVars.nBins_dist_z);

						// Storing probe positions in x2, y2, z2 variables for easy calculations
						// and checking the distance between the probe and the atoms
						x2 = probePosition.x; y2 = probePosition.y; z2 = probePosition.z;
						distance = sqrt (pow ((x2 - x1), 2) + pow ((y2 - y1), 2) + pow ((z2 - z1), 2));

						// Loop through the distribution struct and check if the calculated distance falls within the range
						for (int o = 0; o < freeVolumeVars.nBins_dist; ++o)
						{
							if (distance > (*freeVolumeDist)[o].binStart_dist && distance <= (*freeVolumeDist)[o].binEnd_dist)
							{
								if (isOccupied[index1d] == 0)
								{
									(*freeVolumeDist)[o].nUnoccupied++;
								}
								else if (isOccupied[index1d] == 1)
								{
									(*freeVolumeDist)[o].nOccupied++;
								}
							}
						}
						probePosition.z += freeVolumeVars.delDistance;
					}
					probePosition.y += freeVolumeVars.delDistance;
				}
				probePosition.x += freeVolumeVars.delDistance;
			}
		}
	}

	printf("\n");
}

void initializeFreeVolumeDistribution (FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars)
{
	(*freeVolumeDist)[0].binStart_dist = freeVolumeVars.binStart_dist;
	(*freeVolumeDist)[0].binEnd_dist = (*freeVolumeDist)[0].binStart_dist + freeVolumeVars.binSize_dist;
	(*freeVolumeDist)[0].nUnoccupied = 0;
	(*freeVolumeDist)[0].nOccupied = 0;

	// Setting bin values for free volume distribution
	for (int i = 1; i < freeVolumeVars.nBins_dist; ++i)
	{
		(*freeVolumeDist)[i].binStart_dist = (*freeVolumeDist)[i - 1].binEnd_dist;
		(*freeVolumeDist)[i].binEnd_dist = (*freeVolumeDist)[i].binStart_dist + freeVolumeVars.binSize_dist;
		(*freeVolumeDist)[i].nUnoccupied = 0;
		(*freeVolumeDist)[i].nOccupied = 0;
	}
}

void initializeNBins (FREEVOLUME_VARS *freeVolumeVars)
{
	(*freeVolumeVars).delDistance = 0.25 * (*freeVolumeVars).currentProbeSize;
	(*freeVolumeVars).nBins_dist_x = (int) ((*freeVolumeVars).xLength / (*freeVolumeVars).delDistance);
	(*freeVolumeVars).nBins_dist_y = (int) ((*freeVolumeVars).yLength / (*freeVolumeVars).delDistance);
	(*freeVolumeVars).nBins_dist_z = (int) ((*freeVolumeVars).zLength / (*freeVolumeVars).delDistance);
}

void resetFreeVolumeDistCounts (FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars)
{
	for (int k = 0; k < freeVolumeVars.nBins_dist; ++k)
	{
		(*freeVolumeDist)[k].nUnoccupied = 0;
		(*freeVolumeDist)[k].nOccupied = 0;
	}	
}

void printFreeVolumeDistribution (FREEVOLUME_DISTRIBUTION *freeVolumeDist, int atomType, float probeSize, int nBins)
{
	char *freeVolumeLogfilename;
	freeVolumeLogfilename = (char *) malloc (50 * sizeof (char));
	FILE *freeVolumeLogfile;
	sprintf (freeVolumeLogfilename, "logs/freeVolume_%d_%.5f.log", atomType, probeSize);
	freeVolumeLogfile = fopen (freeVolumeLogfilename, "a");

	for (int i = 0; i < nBins; ++i)
	{
		fprintf(freeVolumeLogfile, "%f %f %d %d %f\n", freeVolumeDist[i].binStart_dist, freeVolumeDist[i].binEnd_dist, freeVolumeDist[i].nOccupied, freeVolumeDist[i].nUnoccupied, ((float)freeVolumeDist[i].nUnoccupied / ((float)freeVolumeDist[i].nUnoccupied + (float)freeVolumeDist[i].nOccupied)));
	}
	free (freeVolumeLogfilename);
	fclose (freeVolumeLogfile);
}

void computeFreeVolume (FREEVOLUME_VARS freeVolumeVars, DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int nThreads)
{
	DATA_ATOMS probePosition;
	freeVolumeVars.currentProbeSize = freeVolumeVars.minProbeSize;

	FREEVOLUME_DISTRIBUTION *freeVolumeDist;
	freeVolumeDist = (FREEVOLUME_DISTRIBUTION *) malloc (freeVolumeVars.nBins_dist * sizeof (FREEVOLUME_DISTRIBUTION));

	initializeFreeVolumeDistribution (&freeVolumeDist, freeVolumeVars);

	float distance, x1, y1, z1, x2, y2, z2;
	long long int index1d, arraySize;

	omp_set_num_threads (nThreads);

	// Probe size is set. It'll vary in a loop, from minimum to maximum size
	for (int i = 0; i < freeVolumeVars.nBins_probeSweep; ++i)
	{
		initializeNBins (&freeVolumeVars);

		probePosition.x = dumpfile.xlo; probePosition.y = dumpfile.ylo; probePosition.z = dumpfile.zlo;

		int *isOccupied;
		arraySize = (long long int)freeVolumeVars.nBins_dist_x * (long long int)freeVolumeVars.nBins_dist_y * (long long int)freeVolumeVars.nBins_dist_z;
		isOccupied = (int *) malloc (arraySize * sizeof (int));
		isOccupied = computeFreeVolume_checkOccupation (i, dumpfile, dumpAtoms, freeVolumeVars, vwdSize, nThreads);

		// Looping through atom types
		// Getting radial distribution of free volume from every atom type
		for (int j = 0; j < entries.nLines_freeVolumeconfig; ++j)
		{
			computeFreeVolume_getDistribution (i, j, &freeVolumeDist, freeVolumeVars, isOccupied, dumpfile, dumpAtoms, freeVolumeconfig, nThreads);
			printFreeVolumeDistribution (freeVolumeDist, freeVolumeconfig[j].atom1, freeVolumeVars.currentProbeSize, freeVolumeVars.nBins_dist);
			resetFreeVolumeDistCounts (&freeVolumeDist, freeVolumeVars);
		}

		freeVolumeVars.currentProbeSize += freeVolumeVars.delProbeSize;
		if (freeVolumeVars.currentProbeSize > freeVolumeVars.maxProbeSize)
			break;

		free (isOccupied);
	}

	free (freeVolumeDist);
}

typedef struct bounds
{
	float lo, hi;
} BOUNDS;

DATA_ATOMS *assignPeaks (DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, BOUNDS *peakInfo, int nPeaks, CONFIG *inputVectors, int nThreads)
{
	DATA_ATOMS *dumpAtomsMod;
	dumpAtomsMod = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));

	float xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo), simVolume = (xDist * yDist * zDist), bondDensity;
	float xDistHalf = (xDist / 2), yDistHalf = (yDist / 2), zDistHalf = (zDist / 2);
	float x_translated, y_translated, z_translated;
	float *lowerBounds, *upperBounds;

	float x1, y1, z1, x2, y2, z2, distance;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		dumpAtomsMod[i].id = dumpAtoms[i].id; dumpAtomsMod[i].atomType = dumpAtoms[i].atomType; dumpAtomsMod[i].molType = 0; dumpAtomsMod[i].x = dumpAtoms[i].x; dumpAtomsMod[i].y = dumpAtoms[i].y; dumpAtomsMod[i].z = dumpAtoms[i].z;
	}

	omp_set_num_threads (nThreads);
	#pragma omp parallel for
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom2Type == inputVectors[0].atom1 && bonds[i].atom1Type == inputVectors[0].atom2))
		{
			bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType; bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;
			bonds[i].x1 = dumpAtoms[bonds[i].atom1 - 1].x; bonds[i].y1 = dumpAtoms[bonds[i].atom1 - 1].y; bonds[i].z1 = dumpAtoms[bonds[i].atom1 - 1].z;
			bonds[i].x2 = dumpAtoms[bonds[i].atom2 - 1].x; bonds[i].y2 = dumpAtoms[bonds[i].atom2 - 1].y; bonds[i].z2 = dumpAtoms[bonds[i].atom2 - 1].z;

			bonds[i].xc = findBondCenter (bonds[i].x1, dumpAtoms[bonds[i].atom1 - 1].ix, bonds[i].x2, dumpAtoms[bonds[i].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
			bonds[i].yc = findBondCenter (bonds[i].y1, dumpAtoms[bonds[i].atom1 - 1].iy, bonds[i].y2, dumpAtoms[bonds[i].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
			bonds[i].zc = findBondCenter (bonds[i].z1, dumpAtoms[bonds[i].atom1 - 1].iz, bonds[i].z2, dumpAtoms[bonds[i].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

			for (int j = 0; j < datafile.nBonds; ++j)
			{
				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom2Type == inputVectors[1].atom1 && bonds[j].atom1Type == inputVectors[1].atom2))
				{
					bonds[j].atom1Type = (int) dumpAtoms[bonds[j].atom1 - 1].atomType; bonds[j].atom2Type = (int) dumpAtoms[bonds[j].atom2 - 1].atomType;
					bonds[j].x1 = dumpAtoms[bonds[j].atom1 - 1].x; bonds[j].y1 = dumpAtoms[bonds[j].atom1 - 1].y; bonds[j].z1 = dumpAtoms[bonds[j].atom1 - 1].z;
					bonds[j].x2 = dumpAtoms[bonds[j].atom2 - 1].x; bonds[j].y2 = dumpAtoms[bonds[j].atom2 - 1].y; bonds[j].z2 = dumpAtoms[bonds[j].atom2 - 1].z;

					bonds[j].xc = findBondCenter (bonds[j].x1, dumpAtoms[bonds[j].atom1 - 1].ix, bonds[j].x2, dumpAtoms[bonds[j].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
					bonds[j].yc = findBondCenter (bonds[j].y1, dumpAtoms[bonds[j].atom1 - 1].iy, bonds[j].y2, dumpAtoms[bonds[j].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
					bonds[j].zc = findBondCenter (bonds[j].z1, dumpAtoms[bonds[j].atom1 - 1].iz, bonds[j].z2, dumpAtoms[bonds[j].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

					x_translated = translatePeriodicDistance (bonds[i].xc, bonds[j].xc, xDistHalf);
					y_translated = translatePeriodicDistance (bonds[i].yc, bonds[j].yc, yDistHalf);
					z_translated = translatePeriodicDistance (bonds[i].zc, bonds[j].zc, zDistHalf);

					distance = sqrt (pow ((bonds[i].xc - x_translated), 2) + pow ((bonds[i].yc - y_translated), 2) + pow ((bonds[i].zc - z_translated), 2));

					for (int k = 0; k < nPeaks; ++k)
					{
						if ((distance <= peakInfo[k].hi) && (distance > peakInfo[k].lo))
						{
							if ((dumpAtomsMod[bonds[j].atom1 - 1].molType == 0) || (dumpAtomsMod[bonds[j].atom1 - 1].molType > (k + 1)))
								dumpAtomsMod[bonds[j].atom1 - 1].molType = (k + 1);

							if ((dumpAtomsMod[bonds[j].atom2 - 1].molType == 0) || (dumpAtomsMod[bonds[j].atom2 - 1].molType > (k + 1)))
								dumpAtomsMod[bonds[j].atom2 - 1].molType = (k + 1);
						}
					}
				}
			}
		}
	}

	return dumpAtomsMod;
}

typedef struct HBondNetwork
{
	int bondPresent, bondAbsent;
} HBONDNETWORK;

void analyzeHBondNetwork (DATA_ATOMS *dumpAtomsMod, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, CONFIG *inputVectors, BOUNDS *peakInfo, int nPeaks, float peakHBondPosition, int currentDumpstep, int nThreads)
{
	float x_translated, y_translated, z_translated, distance;
	float xDistHalf = (dumpfile.xhi - dumpfile.xlo) / 2, yDistHalf = (dumpfile.yhi - dumpfile.ylo) / 2, zDistHalf = (dumpfile.zhi - dumpfile.zlo) / 2;

	// Initializing variables to store bond network information
	int *nHBonds, *nMolecules;
	nHBonds = (int *) calloc ((nPeaks - 1), sizeof (int));
	nMolecules = (int *) calloc (nPeaks + 1, sizeof (int));

	// Writing the output to logfile
	FILE *hBondNetwork_logfile;
	hBondNetwork_logfile = fopen ("hBondNetwork_logs/hBondNetwork.log", "a");

	int loopCounter = 0, index;

	char *connectivityInfo_filename;
	connectivityInfo_filename = (char *) malloc (200 * sizeof (char));

	// Iterating through all atom pairs to check the H-O distance
	omp_set_num_threads (nThreads); printf("\n");
	#pragma omp parallel for
	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		loopCounter++;
		if ((loopCounter % 1000) == 0)
		{
			printf("Checking h-bond networks... %2.3f %%\r", (float) loopCounter * 100.0 / (float) datafile.nAtoms);
			fflush (stdout);
		}

		if (dumpAtomsMod[i].atomType == inputVectors[1].atom1 || dumpAtomsMod[i].atomType == inputVectors[1].atom2)
		{
			for (int j = 0; j < datafile.nAtoms; ++j)
			{
				if (dumpAtomsMod[j].atomType == inputVectors[1].atom1 || dumpAtomsMod[j].atomType == inputVectors[1].atom2)
				{
					x_translated = translatePeriodicDistance (dumpAtomsMod[i].x, dumpAtomsMod[j].x, xDistHalf);
					y_translated = translatePeriodicDistance (dumpAtomsMod[i].y, dumpAtomsMod[j].y, yDistHalf);
					z_translated = translatePeriodicDistance (dumpAtomsMod[i].z, dumpAtomsMod[j].z, zDistHalf);

					distance = sqrt (pow ((dumpAtomsMod[i].x - x_translated), 2) + pow ((dumpAtomsMod[i].y - y_translated), 2) + pow ((dumpAtomsMod[i].z - z_translated), 2));

					// Counting the number of H-bonds between the layers
					if ((distance < peakHBondPosition) && (abs (dumpAtomsMod[i].molType - dumpAtomsMod[j].molType) == 1))
					{
						// If molType of [i] == 2 && molType of [j] == 1,
						// then the pair belongs to first and second layer
						// (i.e), if (molType i - molType j) == 1, then store the bondPresent and bondAbsent stats to [molType i - 2] array index

						if (dumpAtomsMod[i].molType > dumpAtomsMod[j].molType)
						{
							index = dumpAtomsMod[i].molType - 2;

							if (index >= 0)
							{
								snprintf (connectivityInfo_filename, 200, "hBondNetwork_logs/%d-%d-%d_hBond.log", dumpAtomsMod[j].id, dumpAtomsMod[i].id, index);
								FILE *connectivityInfo_file;
								connectivityInfo_file = fopen (connectivityInfo_filename, "a");
								fprintf(connectivityInfo_file, "%d\n", currentDumpstep);
								fclose (connectivityInfo_file);
							}
						}
						else
						{
							index = dumpAtomsMod[j].molType - 2;

							if (index >= 0)
							{
								snprintf (connectivityInfo_filename, 200, "hBondNetwork_logs/%d-%d-%d_hBond.log", dumpAtomsMod[i].id, dumpAtomsMod[j].id, index);
								FILE *connectivityInfo_file;
								connectivityInfo_file = fopen (connectivityInfo_filename, "a");
								fprintf(connectivityInfo_file, "%d\n", currentDumpstep);
								fclose (connectivityInfo_file);
							}
						}

						nHBonds[index]++;
					}

					// creating logfile containing the currentDumpstep and connectivity information
					// connectivity information: 1 => connected; 0 => not connected.
					// it is not necessary to print '0' to the log file
					// if ((distance < peakHBondPosition))
					// {
					// 	if (dumpAtomsMod[i].id > dumpAtomsMod[j].id)
					// 	{
					// 		snprintf (connectivityInfo_filename, 200, "hBondNetwork_logs/%d-%d_hBond.log", dumpAtomsMod[j].id, dumpAtomsMod[i].id);
					// 		FILE *connectivityInfo_file;
					// 		connectivityInfo_file = fopen (connectivityInfo_filename, "a");
					// 		fprintf(connectivityInfo_file, "%d\n", currentDumpstep);
					// 		fclose (connectivityInfo_file);
					// 	}
					// 	else
					// 	{
					// 		snprintf (connectivityInfo_filename, 200, "hBondNetwork_logs/%d-%d_hBond.log", dumpAtomsMod[i].id, dumpAtomsMod[j].id);
					// 		FILE *connectivityInfo_file;
					// 		connectivityInfo_file = fopen (connectivityInfo_filename, "a");
					// 		fprintf(connectivityInfo_file, "%d\n", currentDumpstep);
					// 		fclose (connectivityInfo_file);
					// 	}
					// }
				}
			}

			// Counting the number of molecules in every shell
			nMolecules[ dumpAtomsMod[i].molType ]++;
		}
	}

	// In output, the final result is divided by 2.
	// Because every bond was counted twice.
	// (i.e) once when counting first to second layer,
	// another time when counting from second to first layer.
	for (int i = 0; i < (nPeaks - 1); ++i)
		fprintf(hBondNetwork_logfile, "%.0f, ", (float) (nHBonds[i] / 2));

	for (int i = 0; i < nPeaks; ++i)
		fprintf(hBondNetwork_logfile, "%d, ", nMolecules[i]);

	fprintf(hBondNetwork_logfile, "%d\n", nMolecules[nPeaks]);

	fclose (hBondNetwork_logfile);
}

float *findHBondCorrelation (const char *filename, int nTimeframes, int nThreads)
{
	char lineString[1000], *hBondLogfilename;
	hBondLogfilename = (char *) malloc (100 * sizeof (char));
	snprintf (hBondLogfilename, 100, "./hBondNetwork_logs/%s", filename);

	FILE *hBondLogfile;
	hBondLogfile = fopen (hBondLogfilename, "r");

	// *inputValue contains information about the presence/absence of h-bond
	int lineInt;
	int mean = 0, covariance = 0;

	float *hBondInfo;
	hBondInfo = (float *) calloc (nTimeframes, sizeof (float));

	// Correlation variable
	float *correlation;
	correlation = (float *) calloc (nTimeframes, sizeof (float));

	// A value of '1' is assigned when h-bond is present
	while ((fgets (lineString, 1000, hBondLogfile) != NULL))
	{
		sscanf (lineString, "%d", &lineInt);
		hBondInfo[lineInt] = 1.0;
	}

	// Calculate H-bond correlation (not the classic autocorrelation function)
	// 'i' represents timelag
	for (int i = 0; i < nTimeframes; ++i)
	{
		// Iterating through all the elements
		for (int j = 0; j < (nTimeframes - i); ++j)
		{
			// Calculate h(0).h(t) and find the average
			correlation[i] += hBondInfo[j] * hBondInfo[j + i];
		}
	}

	for (int i = 0; i < nTimeframes; ++i)
	{
		correlation[i] /= nTimeframes;
	}

	fclose (hBondLogfile);
	return correlation;
}

int countNTimeframes (FILE *inputDumpFile)
{
	rewind (inputDumpFile);

	int nLines = 0, nTimeframes;
	char lineString[1000];

	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);

	while ((fgets (lineString, 1000, inputDumpFile) != NULL))
	{
		nLines++;
	}

	rewind (inputDumpFile);

	printf("Number of lines in the input dumpfile: %d\n", nLines);
	nTimeframes = nLines / (dumpfile.nAtoms + 9);
	printf("Number of timeframes: %d\n", nTimeframes);

	return nTimeframes;
}

void calculateSumCorrelation (float *correlation, float **sumCorrelation, int nTimeframes)
{
	for (int i = 0; i < nTimeframes; ++i)
	{
		(*sumCorrelation)[i] += correlation[i];
	}
}

float *calculateAverageCorrelation (float *sumCorrelation, int nFiles, int nTimeframes)
{
	float *avgCorrelation;
	avgCorrelation = (float *) malloc (nTimeframes * sizeof (float));

	for (int i = 0; i < nTimeframes; ++i)
	{
		avgCorrelation[i] = sumCorrelation[i] / nFiles;
	}

	return avgCorrelation;
}

void computeHBondCorrelation2 (const char *fileTemplate, int nTimeframes, int nThreads)
{
	// Create two float arrays, with size equal to the number of timesteps in dump file
	float *avgCorrelation, *sumCorrelation, *correlation;
	sumCorrelation = (float *) calloc (nTimeframes, sizeof (float));
	correlation = (float *) malloc (nTimeframes * sizeof (float));

	struct dirent *filePointer;
	DIR *parentDirectory;
	parentDirectory = opendir ("./hBondNetwork_logs/");

	int nFiles = 0;

	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileTemplate))
		{
			printf("%s\n", filePointer -> d_name);
			correlation = findHBondCorrelation (filePointer -> d_name, nTimeframes, nThreads);
			calculateSumCorrelation (correlation, &sumCorrelation, nTimeframes);
			nFiles++;
		}
	}

	avgCorrelation = calculateAverageCorrelation (sumCorrelation, nFiles, nTimeframes);

	char *bondCorrelation_Filename;
	bondCorrelation_Filename = (char *) malloc (100 * sizeof (char));
	snprintf (bondCorrelation_Filename, 100, "HBond%s.correlation", fileTemplate);

	FILE *bondCorrelation_File;
	bondCorrelation_File = fopen (bondCorrelation_Filename, "w");

	for (int i = 0; i < nTimeframes; ++i)
	{
		fprintf(bondCorrelation_File, "%f\n", log (avgCorrelation[i] / avgCorrelation[0]));
	}

	fclose (bondCorrelation_File);
}

void computeHBondCorrelation (FILE *inputDumpFile, int nThreads)
{
	// Check the number of timesteps in dump file
	int nTimeframes = countNTimeframes (inputDumpFile);

	// Calculate correlation between 1st and 2nd layers (0_hBond.log)
	computeHBondCorrelation2 ("0_hBond.log", nTimeframes, nThreads);
	// Calculate correlation between 2nd and 3rd layers (1_hBond.log)
	computeHBondCorrelation2 ("1_hBond.log", nTimeframes, nThreads);
	// Calculate correlation between 3rd and 4th layers (2_hBond.log)
	computeHBondCorrelation2 ("2_hBond.log", nTimeframes, nThreads);
}

void computeHBonding (DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, BOUNDS *peakInfo, int nPeaks, CONFIG *inputVectors, NLINES_CONFIG entries, float peakHBondPosition, int currentDumpstep, int nThreads)
{
	DATA_ATOMS *dumpAtomsMod;
	dumpAtomsMod = assignPeaks (dumpAtoms, bonds, datafile, dumpfile, peakInfo, nPeaks, inputVectors, nThreads);

	// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// // Checking the previous assignment
	// // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// int nZeroth = 0, nFirst = 0, nSecond = 0, nThird = 0, nFourth = 0;
	// for (int i = 0; i < datafile.nAtoms; ++i)
	// {
	// 	// printf("%d %d %d\n", dumpAtomsMod[i].id, dumpAtomsMod[i].atomType, dumpAtomsMod[i].molType); usleep (10000);
	// 	if (dumpAtomsMod[i].molType == 0)
	// 		nZeroth++;
	// 	if (dumpAtomsMod[i].molType == 1)
	// 		nFirst++;
	// 	if (dumpAtomsMod[i].molType == 2)
	// 		nSecond++;
	// 	if (dumpAtomsMod[i].molType == 3)
	// 		nThird++;
	// 	if (dumpAtomsMod[i].molType == 4)
	// 		nFourth++;
	// }
	// printf("Zeroth: %d\nFirst: %d\nSecond: %d\nThird: %d\nFourth: %d\n", nZeroth, nFirst, nSecond, nThird, nFourth);

	analyzeHBondNetwork (dumpAtomsMod, datafile, dumpfile, inputVectors, peakInfo, nPeaks, peakHBondPosition, currentDumpstep, nThreads);
}

BOUNDS *getHBondPeakInformation (int *nPeaks)
{
	printf("How many peaks to assign?: "); scanf ("%d", &(*nPeaks)); printf("\n");
	BOUNDS *peakInfo;
	peakInfo = (BOUNDS *) malloc ((*nPeaks) * sizeof (BOUNDS));

	for (int i = 0; i < (*nPeaks); ++i)
	{
		printf("Enter the lower bounds for peak %d: ", i + 1); scanf ("%f", &peakInfo[i].lo); printf("\n");
		printf("Enter the upper bounds for peak %d: ", i + 1); scanf ("%f", &peakInfo[i].hi); printf("\n");
	}

	return peakInfo;
}

float getHBondPeakPosition ()
{
	float thresholdHBondDistance;
	printf("Enter the threshold distance to consider for H-bonding: "); scanf ("%f", &thresholdHBondDistance);
	return thresholdHBondDistance;
}

void initializeHBondNetworkLogfile ()
{
	system ("rm hBondNetwork_logs/*");
	FILE *hBondNetwork_logfile;
	hBondNetwork_logfile = fopen ("hBondNetwork_logs/hBondNetwork.log", "w");

	// Writing the first line
	fprintf(hBondNetwork_logfile, "N firstToSecond, N secondToThird, N thirdToFourth, N Mols Zeroth, N Mols First, N Mols Second, N Mols Third, N Mols Fourth\n");

	fclose (hBondNetwork_logfile);
}

void processLAMMPSTraj (FILE *inputDumpFile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int nThreads)
{
	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDumpFile);

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Variables for plotting the distribution
	DIST_VAR plotVars;

	// Calculating the maximum distance between the two farthest points in the simulation box
	float hyp1, xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo);
	hyp1 = sqrt ((xDist * xDist) + (zDist * zDist));
	plotVars.maxDist = sqrt ((hyp1 * hyp1) + (yDist * yDist));

	// Setting the number of bins across distance, degree, and OOP; based on the set bin size. These bin sizes can be adjusted for a smoother distribution curve
	plotVars.binSize_dist = 1; plotVars.binSize_OOP = 0.01; plotVars.binSize_deg = 3;
	plotVars.nBins_dist = (((int) plotVars.maxDist) / (int) plotVars.binSize_dist) + 1; plotVars.nBins_OOP = (int) ((1 + 0.5) / plotVars.binSize_OOP) + 1; plotVars.nBins_deg = (180 / (int) plotVars.binSize_deg) + 1;

	// [degrees][distance] and [oop][distance]
	plotVars.size_degrees = plotVars.nBins_dist * plotVars.nBins_deg; plotVars.size_oop = plotVars.nBins_dist * plotVars.nBins_OOP;
	DISTRIBUTION *distribution_OOP, *distribution_degrees;
	distribution_OOP = (DISTRIBUTION *) malloc (plotVars.size_oop * sizeof (DISTRIBUTION));
	distribution_degrees = (DISTRIBUTION *) malloc (plotVars.size_degrees * sizeof (DISTRIBUTION));

	setDistributionZero (&distribution_degrees, plotVars.size_degrees);
	setDistributionZero (&distribution_OOP, plotVars.size_oop);

	// Setting the bounds of bins (dist, OOP, degrees)
	plotVars.binStart_dist = 0;
	plotVars.binStart_OOP = -0.5;
	plotVars.binStart_deg = 0;
	plotVars.binEnd_dist = plotVars.maxDist;
	plotVars.binEnd_OOP = 1.0;
	plotVars.binEnd_deg = 180;
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	char lineString[1000];
	int isTimestep = 0, currentTimestep, currentLine = 0, isFirstTimestep = 1, currentDumpstep = 0;

	unsigned int nElements = 0;
	int isNElementsSet = 0;

	// Datafile struct is used to store dump atom information
	DATA_ATOMS *dumpAtoms, *dumpAtoms_translated;
	dumpAtoms = (DATA_ATOMS *) malloc (dumpfile.nAtoms * sizeof (DATA_ATOMS));
	dumpAtoms_translated = (DATA_ATOMS *) malloc (dumpfile.nAtoms * sizeof (DATA_ATOMS));

	ORDERPARAMETER *allData_array;

	printf("\n");

	// bondRDF variable
	int RDFcounter = 0;
	float *bondRDF, binSize_dist_RDF = 0.3;
	int nBins_dist_RDF = (int) (plotVars.maxDist / binSize_dist_RDF);
	bondRDF = (float *) malloc (nBins_dist_RDF * sizeof (float));

	// Free volume calculations - variable set
	FREEVOLUME_VARS freeVolumeVars;

	// Assigning peaks
	BOUNDS *peakInfo; int nPeaks; float peakHBondPosition;

	// Reading and processing dump information
	while (fgets (lineString, 1000, inputDumpFile) != NULL)
	{
		// Executed at the end of first timeframe
		if (currentLine == (9 + dumpfile.nAtoms) && nElements == 0)
		{
			nElements = getNElements (datafile, dumpAtoms, bonds, inputVectors);
			plotVars.nElements = nElements;
			printf("Allocating memory for %lu elements...\n", nElements);
			allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));
			printf("Memory allocated successfully...\n");

			freeVolumeVars = getFreeVolumeVars (dumpfile);

			// Gathering peak information for analysing H-bonds
			peakInfo = getHBondPeakInformation (&nPeaks);
			peakHBondPosition = getHBondPeakPosition ();
			initializeHBondNetworkLogfile ();
			printf("HBond network file initialized successfully...\n");
			fflush (stdout);
		}

		// Main processing loop
		if ((currentDumpstep > 2) && (nElements > 0) && (currentLine == 2))
		{
			sscanf (lineString, "%d", &currentTimestep);
			printf("Scanning timestep: %d...               \n", currentTimestep);
			fflush (stdout); 

			// computeBondRDF (dumpAtoms, datafile, dumpfile, bonds, inputVectors, plotVars, nThreads, binSize_dist_RDF, &bondRDF, &RDFcounter, currentTimestep);

			// Checking bond orientation
			// allData_array = computeOrderParameter (dumpAtoms, dumpfile, datafile, bonds, inputVectors, currentTimestep, nElements);
			// computeDistribution_OOP (allData_array, plotVars, &distribution_OOP, nThreads);
			// computeDistribution_theta (allData_array, plotVars, &distribution_degrees, nThreads);

			// Calculating free volume distribution
			// computeFreeVolume (freeVolumeVars, dumpAtoms, dumpfile, freeVolumeconfig, vwdSize, entries, nThreads);

			// Computing H bond lifetime and frequency
			computeHBonding (dumpAtoms, bonds, datafile, dumpfile, peakInfo, nPeaks, inputVectors, entries, peakHBondPosition, currentDumpstep, nThreads);

			isTimestep = 0;
		}

		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			isTimestep = 1;
			currentDumpstep++;
			currentLine = 1;
			isNElementsSet = 0;
		}

		if (currentLine > 9 && currentLine < (9 + dumpfile.nAtoms))
		{	
			sscanf (lineString, "%d %d %f %f %f %*f %*f %*f %d %d %d\n",
				&dumpAtoms[currentLine - 10].id,
				&dumpAtoms[currentLine - 10].atomType,
				&dumpAtoms[currentLine - 10].x,
				&dumpAtoms[currentLine - 10].y,
				&dumpAtoms[currentLine - 10].z,
				&dumpAtoms[currentLine - 10].ix,
				&dumpAtoms[currentLine - 10].iy,
				&dumpAtoms[currentLine - 10].iz);
		}

		currentLine++;
	}

	printDistribution_OOP (distribution_OOP, plotVars);
	printDistribution_degrees (distribution_degrees, plotVars);
	printBondRDF (bondRDF, RDFcounter, nBins_dist_RDF, binSize_dist_RDF);
}

int main(int argc, char const *argv[])
{
	system ("mkdir logs");
	system ("mkdir bondRDF_logs");
	system ("mkdir hBondNetwork_logs");

	long number_of_processors = sysconf(_SC_NPROCESSORS_ONLN);
	int nThreads = (int) number_of_processors - 1;

	FILE *inputDumpFile, *inputDataFile, *inputConfigFile, *inputFreevolumeConfigFile, *inputVDWConfigFile, *inputHBondConfigFile;
	char *inputDumpFilename, *inputDataFilename, *inputConfigFilename, *inputFreevolumeConfigFilename, *inputVDWConfigFilename, *inputHBondConfigFilename;

	printf("%s\n", "Looking for LAMMPS trajectory file...");
	inputDumpFilename = getInputFileName ();
	printf("%s\n", "Looking for LAMMPS data file...");
	inputDataFilename = getInputFileName ();
	printf("%s\n", "Looking for input config file (for OOP/bondRDF calculations)...");
	inputConfigFilename = getInputFileName ();
	printf("%s\n", "Looking for config file for free volume calculations...");
	inputFreevolumeConfigFilename = getInputFileName ();
	printf("%s\n", "Looking for config file containing VWD radii...");
	inputVDWConfigFilename = getInputFileName ();
	printf("%s\n", "Looking for config file to compute H bonding...");
	inputHBondConfigFilename = getInputFileName ();

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

	// processLAMMPSTraj (inputDumpFile, datafile, bonds, inputVectors, freeVolumeconfig, vwdSize, entries, nThreads);

	// Checking h-bond lifetime correlation function from the saved logfiles
	computeHBondCorrelation (inputDumpFile, nThreads);

	return 0;
}