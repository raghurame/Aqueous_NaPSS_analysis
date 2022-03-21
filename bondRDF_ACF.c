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
#include "fileHandling.h"
#include "generalUtilities.h"
#include "bondRDF_ACF.h"

/*
 * Arguments to pass:
 * ~~~~~~~~~~~~~~~~~
 * argv[0] = ./program
 * argv[1] = main LAMMPS traj dump filename
 * argv[2] = threshold distance to consider from the first timeframe
 *
 */

long long int findNlines (char *logFilename)
{
	long long int nLinesTotal = 0;
	char lineString[1000];
	FILE *inputFile;
	inputFile = fopen (logFilename, "r");

	while ((fgets (lineString, 1000, inputFile) != NULL))
	{
		nLinesTotal++;
	}

	fclose (inputFile);
	return nLinesTotal;
}

void storeLogfileInfo (char *logFilename, ALL_DATA_BONDRDF_ACF **fullData, int nLines, int currentTrajCount)
{
	FILE *logFile;
	logFile = fopen (logFilename, "r");

	long long int index1d = 0;
	int currentLine = 0;
	char lineString[1000];

	while ((fgets (lineString, 1000, logFile) != NULL))
	{
		// currentTrajCount is the row
		// currentLine is the column
		// The whole data from input traj is stored in a single row. 
		index1d = getIndex1d (currentTrajCount, currentLine, nLines);
		sscanf (lineString, "%f %f %f %f %f %f %f %f %f %f\n", &(*fullData)[index1d].x1, &(*fullData)[index1d].y1, &(*fullData)[index1d].z1, &(*fullData)[index1d].x2, &(*fullData)[index1d].y2, &(*fullData)[index1d].z2, &(*fullData)[index1d].x3, &(*fullData)[index1d].y3, &(*fullData)[index1d].z3, &(*fullData)[index1d].distance);
		currentLine++;
	}
}

LOGFILES_VARIABLES_BONDRDF_ACF openLogFiles (FILE *mainDumpfile, ALL_DATA_BONDRDF_ACF **fullData)
{
	char lineString[1000], *logFilename;
	int isTimeline = 0, currentTimeframe;

	LOGFILES_VARIABLES_BONDRDF_ACF fileVars;
	fileVars.currentTrajCount = 0;
	fileVars.nLinesTotal = 0;

	logFilename = (char *) malloc (50 * sizeof (char));

	while ((fgets (lineString, 1000, mainDumpfile) != NULL))
	{
		if (isTimeline)
		{
			sscanf (lineString, "%d", &currentTimeframe);
			snprintf (logFilename, 50, "bondRDF_logs/%d.rdf", currentTimeframe);

			if (isFileExists (logFilename))
			{
				fileVars.nLines = findNlines (logFilename);
				fileVars.nLinesTotal += fileVars.nLines;
				(*fullData) = (ALL_DATA_BONDRDF_ACF *) realloc ((*fullData), fileVars.nLinesTotal * sizeof (ALL_DATA_BONDRDF_ACF));
				storeLogfileInfo (logFilename, &(*fullData), fileVars.nLines, fileVars.currentTrajCount);
				printf("reading file: %s     \r", logFilename);
				fflush (stdout);
				fileVars.currentTrajCount++;

				if (fileVars.currentTrajCount >= 1000)
				{
					goto earlyExit;
				}
			}

			isTimeline = 0;
		}

		if (strstr (lineString, "ITEM: TIMESTEP"))
			isTimeline = 1;
	}

	earlyExit:

	return fileVars;
}

// printACF (i, originalDistance, acf, fileVars.currentTrajCount, lowerLimit, upperLimit);
void printACF (int i, float *originalDistance, float *acf, int currentTrajCount, float lowerLimit, float upperLimit)
{
	if ((originalDistance[0] < upperLimit) && (originalDistance[0] > lowerLimit))
	{
		char *acfOutput_filename;
		acfOutput_filename = (char *) malloc (50 * sizeof (char));

		FILE *acfOutput_file;
		snprintf (acfOutput_filename, 50, "bondRDF_processed/processed_%d.rdf", i);
		acfOutput_file = fopen (acfOutput_filename, "w");

		for (int j = 0; j < currentTrajCount; ++j)
		{
			fprintf(acfOutput_file, "%.3f, %.3f\n", originalDistance[j], acf[j]);
		}
		fclose (acfOutput_file);
		free (acfOutput_filename);
	}
}

void computeACF (LOGFILES_VARIABLES_BONDRDF_ACF fileVars, ALL_DATA_BONDRDF_ACF *fullData)
{
	system ("mkdir bondRDF_processed");

	long long int index1d, index1d_2;
	float mean, covariance, covariance_var, *acf, *originalDistance;
	acf = (float *) malloc (fileVars.currentTrajCount * sizeof (float));
	originalDistance = (float *) malloc (fileVars.currentTrajCount * sizeof (float));

	float lowerLimit, upperLimit;

	// 'i' denotes the bond pairs
	for (int i = 0; i < fileVars.nLines; ++i)
	{
		mean = 0; covariance = 0;

		// 'j' denotes the time frame for every bond pair
		// computing mean bond-bond distance
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			mean += fullData[index1d].distance;
		}
		mean /= fileVars.currentTrajCount;

		// computing covariance of bond-bond distance series
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			covariance_var = fullData[index1d].distance - mean;
			covariance += pow ((fullData[index1d].distance - mean), 2);
		}
		covariance /= fileVars.currentTrajCount;

		// subtracting mean from every point
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			originalDistance[j] = fullData[index1d].distance;
			fullData[index1d].distance -= mean;
		}

		// computing autocorrelation
		// in this loop, 'j' is the timelag
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			acf[j] = 0;

			// iterating through all the elements
			for (int k = 0; k < (fileVars.currentTrajCount - j); ++k)
			{
				index1d = getIndex1d (k, i, fileVars.nLines);
				index1d_2 = getIndex1d (k + j, i, fileVars.nLines);
				acf[j] += (fullData[index1d].distance * fullData[index1d_2].distance);
			}

			acf[j] /= fileVars.currentTrajCount;
			acf[j] /= covariance;
			// index1d = getIndex1d (j, i, fileVars.nLines);
		}

		// printing the acf function
		// printing acf function for first peak
		lowerLimit = 0; upperLimit = 5.2;
		printACF (i, originalDistance, acf, fileVars.currentTrajCount, lowerLimit, upperLimit);

		// printing acf function for second peak
		lowerLimit = 5.2; upperLimit = 10;
		printACF (i, originalDistance, acf, fileVars.currentTrajCount, lowerLimit, upperLimit);

		// printing acf function for third peak
		lowerLimit = 10; upperLimit = 15;
		printACF (i, originalDistance, acf, fileVars.currentTrajCount, lowerLimit, upperLimit);

		// printing acf function beyond third peak for comparison
		lowerLimit = 15; upperLimit = 20;
		printACF (i, originalDistance, acf, fileVars.currentTrajCount, lowerLimit, upperLimit);
	}
}

void computeACFOfBondRDF (FILE *inputDumpFile)
{
	ALL_DATA_BONDRDF_ACF *fullData;
	// fullData is allocated for 10 elements
	// 10 elements is an arbitrary number
	// Memory is later reallocated inside openLogFiles () function
	fullData = (ALL_DATA_BONDRDF_ACF *) malloc (10 * sizeof (ALL_DATA_BONDRDF_ACF));

	LOGFILES_VARIABLES_BONDRDF_ACF fileVars;

	fileVars = openLogFiles (inputDumpFile, &fullData);
	computeACF (fileVars, fullData);
}

// int main(int argc, char const *argv[])
// {
// 	long number_of_processors = sysconf(_SC_NPROCESSORS_ONLN);
// 	int nThreads = (int) number_of_processors - 1;

// 	if (argc < 3)
// 	{
// 		printf("Required args:\n~~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = main dump file\nargv[2] = threshold distance to consider\n\n");
// 		exit (1);
// 	}
// 	FILE *mainDumpfile;
// 	mainDumpfile = fopen (argv[1], "r");

// 	float thresholdDistance = atof (argv[2]);
// 	ALL_DATA *fullData;
// 	fullData = (ALL_DATA *) malloc (10 * sizeof (ALL_DATA));

// 	LOGFILES_VARIABLES fileVars;

// 	fileVars = openLogFiles (mainDumpfile, thresholdDistance, &fullData);
// 	computeACF (fileVars, fullData);

// 	return 0;
// }