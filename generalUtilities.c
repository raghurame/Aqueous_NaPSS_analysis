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

float translatePeriodic (float coord, int image, float dimension)
{
	coord += (image * dimension);
	return coord;
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
