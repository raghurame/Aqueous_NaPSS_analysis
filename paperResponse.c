#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <dirent.h>
#include <unistd.h>
#include <stdlib.h>

typedef struct energyDistribution
{
	int count;
	float min, max;
} ENERGY_DISTRIBUTION;

void saveAtomIDs (FILE *inputData, int **atomIDs_O, int *nOs, int ID_O, int **atomIDs_H, int *nHs, int ID_H)
{
	char lineString[2000];
	(*nOs) = 0;
	(*nHs) = 0;

	int isAtomLine = 0, currentType, currentID;

	while (fgets (lineString, 2000, inputData) != NULL)
	{
		if (isAtomLine == 1)
		{
			sscanf (lineString, "%*d %*d %d %*f %*f %*f %*f\n", &currentType);

			if (currentType == ID_O) {
				(*nOs)++; }

			if (currentType == ID_H) {
				(*nHs)++; }
		}

		if (strstr (lineString, "Atoms")) {
			fgets (lineString, 2000, inputData);
			isAtomLine = 1; }

		if (strstr (lineString, "Bonds")) {
			fgets (lineString, 2000, inputData);
			isAtomLine = 0;
			break; }
	}

	rewind (inputData);

	(*atomIDs_O) = (int *) malloc ((*nOs) * sizeof (int));
	(*atomIDs_H) = (int *) malloc ((*nHs) * sizeof (int));

	int counter_O = 0, counter_H = 0;

	while (fgets (lineString, 2000, inputData) != NULL)
	{
		if (isAtomLine == 1)
		{
			sscanf (lineString, "%d %*d %d %*f %*f %*f %*f\n", &currentID, &currentType);

			if (currentType == ID_O) {
				(*atomIDs_O)[counter_O] = currentID;
				counter_O++; }
			else if (currentType == ID_H) {
				(*atomIDs_H)[counter_H] = currentID;
				counter_H++; }
		}

		if (strstr (lineString, "Atoms")) {
			fgets (lineString, 2000, inputData);
			isAtomLine = 1; }

		if (strstr (lineString, "Bonds")) {
			fgets (lineString, 2000, inputData);
			isAtomLine = 0;
			break; }
	}
}

int arrint (int *array, int arraySize, int integer)
{
	for (int i = 0; i < arraySize; ++i)
	{
		if (integer == array[i]) {
			return 1; }
	}

	return 0;
}

void findEnergyBounds (const char *fileName, float *minEnergy, float *maxEnergy, int *atomIDs_O, int nOs, int *atomIDs_H, int nHs)
{
	printf("Computing energy bounds (lower and upper limits)...\n");
	FILE *pairLocalFile;
	char *pipeString, lineString[2000];
	pipeString = (char *) malloc (200 * sizeof (char));

	int countingLine = 0;

	int patom1, patom2;
	float distance, energy;

	if (strstr (fileName, ".xz")) {
		snprintf (pipeString, 200, "xzcat %s", fileName);
		pairLocalFile = popen (pipeString, "r"); }
	else {
		pairLocalFile = fopen (fileName, "r"); }

	// Skipping the first 9 header files
	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, pairLocalFile); }

	// Reading the entries here
	while (fgets (lineString, 2000, pairLocalFile) != NULL)
	{
		// Format in pairLocal file
		// patom1 patom2 distance energy fx fy fz
		sscanf (lineString, "%*d %d %d %f %f %*f %*f %*f\n", &patom1, &patom2, &distance, &energy);

		// Check if patom1 and patom2 belongs to *atomIDs_O or *atomIDs_H
		if ((arrint (atomIDs_O, nOs, patom1) || arrint (atomIDs_O, nOs, patom2)) && (arrint (atomIDs_H, nHs, patom1) || arrint (atomIDs_H, nHs, patom2)))
		{
			if ((*minEnergy) == 0) {
				(*minEnergy) = energy; }

			if ((*maxEnergy) == 0) {
				(*maxEnergy) = energy; }

			if (energy < (*minEnergy)) {
				(*minEnergy) = energy; }

			if (energy > (*maxEnergy)) {
				(*maxEnergy) = energy; }
		}

		countingLine++;

		if ((countingLine % 10000) == 0) {
			printf("Current line...%d                       \r", countingLine);
			fflush (stdout); }
	}
}

void createHistogram (const char *fileName, float minEnergy, float maxEnergy, int *atomIDs_O, int nOs, int *atomIDs_H, int nHs, ENERGY_DISTRIBUTION **engDist, int nBins)
{
	printf("\nCreating histogram...\n");

	int countingLine = 0;

	FILE *input;

	char *pipeString, lineString[2000];
	pipeString = (char *) malloc (200 * sizeof (char));

	int patom1, patom2;
	float distance, energy;

	// Reading the input file
	if (strstr (fileName, ".xz")) {
		snprintf (pipeString, 200, "xzcat %s", fileName);
		input = popen (pipeString, "r"); }
	else {
		input = fopen (fileName, "r"); }

	// Skipping the first 9 header lines in input file
	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, input); }

	// Reading the entries
	while (fgets (lineString, 2000, input) != NULL)
	{
		// Format in pairLocal file
		// patom1 patom2 distance energy fx fy fz
		sscanf (lineString, "%*d %d %d %f %f %*f %*f %*f\n", &patom1, &patom2, &distance, &energy);

		// Check if patom1 and patom2 belongs to *atomIDs_O or *atomIDs_H
		if ((arrint (atomIDs_O, nOs, patom1) || arrint (atomIDs_O, nOs, patom2)) && (arrint (atomIDs_H, nHs, patom1) || arrint (atomIDs_H, nHs, patom2)))
		{
			// Creating an histogram
			for (int i = 0; i < nBins; ++i)
			{
				if (energy < (*engDist)[i].max && energy >= (*engDist)[i].min) {
					(*engDist)[i].count++; }
			}
		}

		countingLine++;

		if ((countingLine % 10000) == 0) {
			printf("Current line...%d                         \r", countingLine);
			fflush (stdout); }
	}
}

ENERGY_DISTRIBUTION *assignBinValues (ENERGY_DISTRIBUTION *engDist, float minEnergy, float maxEnergy, int nBins)
{
	// Allocating bounds for each bins
	float overallEnergyWidth = (maxEnergy - minEnergy), binWidth = overallEnergyWidth / (float) nBins;

	// Allocating 'nBins' points for the distribution
	engDist = (ENERGY_DISTRIBUTION *) malloc (nBins * sizeof (ENERGY_DISTRIBUTION));

	// Assigning min and max values for engDist
	engDist[0].min = minEnergy;
	engDist[0].max = minEnergy + binWidth;
	engDist[0].count = 0;

	for (int i = 1; i < nBins; ++i)
	{
		engDist[i].min = engDist[i - 1].max;
		engDist[i].max = engDist[i].min + binWidth;
		engDist[i].count = 0;
	}

	return engDist;
}

int main(int argc, char const *argv[])
{
	DIR *parentDirectory;
	struct dirent *filePointer;

	parentDirectory = opendir ("./");

	ENERGY_DISTRIBUTION *engDist;

	// Save atom IDs for Os (from PSS) and Hs (from water).
	FILE *inputData;
	inputData = fopen (argv[1], "r");

	int *atomIDs_O, nOs, ID_O = atoi (argv[2]), *atomIDs_H, nHs, ID_H = atoi (argv[3]);

	saveAtomIDs (inputData, &atomIDs_O, &nOs, ID_O, &atomIDs_H, &nHs, ID_H);

	float minEnergy = 0, maxEnergy = 0;
	int nBins = 50;

	int nMaxFilesToCheck = 999999, currentFile = 1;
	while ((filePointer = readdir (parentDirectory)))
	{
		if (strstr (filePointer -> d_name, "pairLocal.") && strstr (filePointer -> d_name, ".dump."))
		{
			printf("\n~~~~~~~~~~~~~~~~\nPROCESSING FILE: %s\n~~~~~~~~~~~~~~~~\n\n", filePointer -> d_name);

			// Find the min and max energy between the pairs
			findEnergyBounds (filePointer -> d_name, &minEnergy, &maxEnergy, atomIDs_O, nOs, atomIDs_H, nHs);
			printf("minEnergy: %f; maxEnergy: %f;\n", minEnergy, maxEnergy);

			currentFile++;
			if (currentFile > nMaxFilesToCheck) {
				goto leaveThisWhileLoop; }
		}
	}

	leaveThisWhileLoop: ;

	// Rewind the directory pointer
	rewinddir (parentDirectory);

	// Assigning bin values
	engDist = assignBinValues (engDist, minEnergy, maxEnergy, nBins);

	currentFile = 1;
	char *individualFileName;
	individualFileName = (char *) malloc (200 * sizeof (char));

	system ("mkdir logs");

	while ((filePointer = readdir (parentDirectory)))
	{
		if (strstr (filePointer -> d_name, "pairLocal.") && strstr (filePointer -> d_name, ".dump."))
		{
			printf("\n~~~~~~~~~~~~~~~~\nPROCESSING FILE: %s\n~~~~~~~~~~~~~~~~\n\n", filePointer -> d_name);

			// Create a histogram and compute the distribution
			createHistogram (filePointer -> d_name, minEnergy, maxEnergy, atomIDs_O, nOs, atomIDs_H, nHs, &engDist, nBins);

			FILE *outputIndividualStats;
			snprintf (individualFileName, 200, "logs/energy.%d.dist", currentFile);
			printf("Printing file %s\n", individualFileName);
			outputIndividualStats = fopen (individualFileName, "w");

			for (int i = 0; i < nBins; ++i)
			{
				fprintf(outputIndividualStats, "%f %f %d\n", engDist[i].min, engDist[i].max, engDist[i].count);
			}

			fclose (outputIndividualStats);

			currentFile++;
			if (currentFile > nMaxFilesToCheck) {
				goto leaveThisWhileLoop2; }
		}
	}

	leaveThisWhileLoop2: ;

	// Printing the energy distribution
	FILE *outputFile;
	outputFile = fopen ("energy.dist", "w");

	for (int i = 0; i < nBins; ++i) {
		fprintf(outputFile, "%f %f %d\n", engDist[i].min, engDist[i].max, engDist[i].count); }

	return 0;
}