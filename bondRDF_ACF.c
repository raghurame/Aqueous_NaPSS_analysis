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

void storeLogfileInfo (char *logFilename, float **mainData, int nLines, int currentTrajCount)
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
		sscanf (lineString, "%*d %*d %*d %*d %f\n", &(*mainData)[index1d]);
		currentLine++;
	}
}

typedef struct openLogFiles_struct
{
	int currentTrajCount, nLines;
	long long int nLinesTotal;
} LOGFILES_VARIABLES;

LOGFILES_VARIABLES openLogFiles (FILE *mainDumpfile, float thresholdDistance, float **mainData)
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
				(*mainData) = (float *) realloc ((*mainData), fileVars.nLinesTotal * sizeof (float));
				storeLogfileInfo (logFilename, &(*mainData), fileVars.nLines, fileVars.currentTrajCount);
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

void computeACF (LOGFILES_VARIABLES fileVars, float *mainData)
{
	long long int index1d, index1d_2;
	float mean, covariance, covariance_var, *acf;
	acf = (float *) malloc (fileVars.currentTrajCount * sizeof (float));
	int 

	// 'i' denotes the bond pairs
	for (int i = 0; i < fileVars.nLines; ++i)
	{
		mean = 0; covariance = 0;

		// 'j' denotes the time frame for every bond pair
		// computing mean bond-bond distance
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			mean += mainData[index1d];
		}
		mean /= fileVars.currentTrajCount;

		// computing covariance of bond-bond distance series
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			covariance_var = mainData[index1d] - mean;
			covariance += pow ((mainData[index1d] - mean), 2);
		}
		covariance /= fileVars.currentTrajCount;

		// computing autocorrelation
		// subtracting mean from every point
		for (int j = 0; j < fileVars.currentTrajCount; ++j)
		{
			index1d = getIndex1d (j, i, fileVars.nLines);
			mainData[index1d] -= mean;
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
				acf[j] += mainData[index1d] * mainData[index1d_2];
			}
			acf[j] /= fileVars.currentTrajCount;
			acf[j] /= covariance;
		}
	}
}

int main(int argc, char const *argv[])
{
	FILE *mainDumpfile;
	mainDumpfile = fopen (argv[1], "r");

	float thresholdDistance = atof (argv[2]);
	float *mainData;
	mainData = (float *) malloc (10 * sizeof (float));

	LOGFILES_VARIABLES fileVars;

	fileVars = openLogFiles (mainDumpfile, thresholdDistance, &mainData);
	computeACF (fileVars, mainData);

	return 0;
}