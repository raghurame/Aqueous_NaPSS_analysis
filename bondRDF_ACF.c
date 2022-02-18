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

/*
 * Arguments to pass:
 * ~~~~~~~~~~~~~~~~~
 * argv[0] = ./program
 * argv[1] = main LAMMPS traj dump filename
 * argv[2] = threshold distance to consider from the first timeframe
 *
 */

int isFileExists (char *inputFilename)
{
	FILE *checking;

	if (checking = fopen (inputFilename, "r"))
	{
		fclose (checking);
		return 1;
	}
	return 0;
}

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

long long int getIndex1d (int i, int j, int width_j)
{
	// width_j is the max value of 'j'
	int index1d = 0;
	index1d = (width_j * i) + j;
	return index1d;
}

typedef struct allData
{
	float x1, y1, z1, x2, y2, z2, x3, y3, z3, distance;
} ALL_DATA;

void storeLogfileInfo (char *logFilename, ALL_DATA **fullData, int nLines, int currentTrajCount)
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

typedef struct openLogFiles_struct
{
	int currentTrajCount, nLines;
	long long int nLinesTotal;
} LOGFILES_VARIABLES;

LOGFILES_VARIABLES openLogFiles (FILE *mainDumpfile, float thresholdDistance, ALL_DATA **fullData)
{
	char lineString[1000], *logFilename;
	int isTimeline = 0, currentTimeframe;

	LOGFILES_VARIABLES fileVars;
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
				(*fullData) = (ALL_DATA *) realloc ((*fullData), fileVars.nLinesTotal * sizeof (ALL_DATA));
				storeLogfileInfo (logFilename, &(*fullData), fileVars.nLines, fileVars.currentTrajCount);
				printf("reading file: %s     \r", logFilename);
				fflush (stdout);
				fileVars.currentTrajCount++;
			}

			isTimeline = 0;
		}

		if (strstr (lineString, "ITEM: TIMESTEP"))
			isTimeline = 1;
	}

	return fileVars;
}

void computeACF (LOGFILES_VARIABLES fileVars, ALL_DATA *fullData)
{
	FILE *sampleOutput;
	sampleOutput = fopen ("acf_sample", "w");

	long long int index1d, index1d_2;
	float mean, covariance, covariance_var, *acf;
	acf = (float *) malloc (fileVars.currentTrajCount * sizeof (float));

	// 'i' denotes the bond pairs
	for (int i = 0; i < fileVars.nLines; ++i)
	{
		mean = 0; covariance = 0;

		// 'j' denotes the time frame for every bond pair
		// computing mean bond-bond distance
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			printf("%f\n", fullData[index1d].distance);
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

		// computing autocorrelation
		// subtracting mean from every point
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			fullData[index1d].distance -= mean;
		}

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
			index1d = getIndex1d (j, i, fileVars.nLines);
			fprintf(sampleOutput, "[%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f], [%.3f, %.3f, %.3f], %.3f, %.3f\n", fullData[index1d].x1, fullData[index1d].y1, fullData[index1d].z1, fullData[index1d].x2, fullData[index1d].y2, fullData[index1d].z2, fullData[index1d].x3, fullData[index1d].y3, fullData[index1d].z3, fullData[index1d].distance, acf[j]);
		}
		printf("\n");
		exit (1);
	}

	fclose (sampleOutput);
}

int main(int argc, char const *argv[])
{
	if (argc < 3)
	{
		printf("Required args:\n~~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = main dump file\nargv[2] = threshold distance to consider\n\n");
		exit (1);
	}
	FILE *mainDumpfile;
	mainDumpfile = fopen (argv[1], "r");

	float thresholdDistance = atof (argv[2]);
	ALL_DATA *fullData;
	fullData = (ALL_DATA *) malloc (10 * sizeof (ALL_DATA));

	LOGFILES_VARIABLES fileVars;

	fileVars = openLogFiles (mainDumpfile, thresholdDistance, &fullData);
	computeACF (fileVars, fullData);

	return 0;
}