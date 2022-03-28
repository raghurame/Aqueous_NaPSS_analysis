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
#include "hBondCorrelation.h"
#include "generalUtilities.h"
#include "fileHandling.h"

DATA_ATOMS *assignPeaks (DATA_ATOMS *dumpAtoms, DATA_ATOMS *dumpAtomsMod, DATA_BONDS *bonds, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, BOUNDS *peakInfo, int nPeaks, CONFIG *inputVectors, int nThreads)
{
	// DATA_ATOMS *dumpAtomsMod;
	// dumpAtomsMod = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));

	float xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo);
	float xDistHalf = (xDist / 2), yDistHalf = (yDist / 2), zDistHalf = (zDist / 2);
	float x_translated, y_translated, z_translated;

	float distance;

	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		dumpAtomsMod[i].id = dumpAtoms[i].id; 
		dumpAtomsMod[i].atomType = dumpAtoms[i].atomType; 
		dumpAtomsMod[i].molType = 0; 
		dumpAtomsMod[i].x = dumpAtoms[i].x; 
		dumpAtomsMod[i].y = dumpAtoms[i].y; 
		dumpAtomsMod[i].z = dumpAtoms[i].z;
	}

	omp_set_num_threads (nThreads);
	#pragma omp parallel for
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom2Type == inputVectors[0].atom1 && bonds[i].atom1Type == inputVectors[0].atom2))
		{
			bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType; 
			bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;

			bonds[i].x1 = dumpAtoms[bonds[i].atom1 - 1].x; 
			bonds[i].y1 = dumpAtoms[bonds[i].atom1 - 1].y; 
			bonds[i].z1 = dumpAtoms[bonds[i].atom1 - 1].z;

			bonds[i].x2 = dumpAtoms[bonds[i].atom2 - 1].x; 
			bonds[i].y2 = dumpAtoms[bonds[i].atom2 - 1].y; 
			bonds[i].z2 = dumpAtoms[bonds[i].atom2 - 1].z;

			bonds[i].xc = findBondCenter (bonds[i].x1, dumpAtoms[bonds[i].atom1 - 1].ix, bonds[i].x2, dumpAtoms[bonds[i].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
			bonds[i].yc = findBondCenter (bonds[i].y1, dumpAtoms[bonds[i].atom1 - 1].iy, bonds[i].y2, dumpAtoms[bonds[i].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
			bonds[i].zc = findBondCenter (bonds[i].z1, dumpAtoms[bonds[i].atom1 - 1].iz, bonds[i].z2, dumpAtoms[bonds[i].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

			for (int j = 0; j < datafile.nBonds; ++j)
			{
				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom2Type == inputVectors[1].atom1 && bonds[j].atom1Type == inputVectors[1].atom2))
				{
					bonds[j].atom1Type = (int) dumpAtoms[bonds[j].atom1 - 1].atomType; 
					bonds[j].atom2Type = (int) dumpAtoms[bonds[j].atom2 - 1].atomType;

					bonds[j].x1 = dumpAtoms[bonds[j].atom1 - 1].x; 
					bonds[j].y1 = dumpAtoms[bonds[j].atom1 - 1].y; 
					bonds[j].z1 = dumpAtoms[bonds[j].atom1 - 1].z;

					bonds[j].x2 = dumpAtoms[bonds[j].atom2 - 1].x; 
					bonds[j].y2 = dumpAtoms[bonds[j].atom2 - 1].y; 
					bonds[j].z2 = dumpAtoms[bonds[j].atom2 - 1].z;

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
	// #pragma omp parallel for
	for (int i = 0; i < datafile.nAtoms; ++i)
	{
		loopCounter++;
		if ((loopCounter % 1000) == 0)
		{
			printf("Checking h-bond networks... %2.3f %%\r", (float) loopCounter * 100.0 / (float) datafile.nAtoms);
			fflush (stdout);
		}

		// The second line in bondRDF.config contains dipoles in solvent
		// This 'if' statement checks if the 'i'th atom corresponds to the solvent atom
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

	free (nHBonds); 
	free (nMolecules);

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

	// for (int i = 0; i < nTimeframes; ++i)
	// {
	// 	correlation[i] /= nTimeframes;
	// }

	fclose (hBondLogfile);
	return correlation;
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
		if (avgCorrelation[i] > 0)
		{
			fprintf(bondCorrelation_File, "%f %f %f %f\n", 
				log (avgCorrelation[i] / avgCorrelation[0]), 
				(avgCorrelation[i] / avgCorrelation[0]), 
				avgCorrelation[i], 
				avgCorrelation[0]);
		}
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
	dumpAtomsMod = (DATA_ATOMS *) malloc (datafile.nAtoms * sizeof (DATA_ATOMS));

	dumpAtomsMod = assignPeaks (dumpAtoms, dumpAtomsMod, bonds, datafile, dumpfile, peakInfo, nPeaks, inputVectors, nThreads);

	analyzeHBondNetwork (dumpAtomsMod, datafile, dumpfile, inputVectors, peakInfo, nPeaks, peakHBondPosition, currentDumpstep, nThreads);

	free (dumpAtomsMod);
}

BOUNDS *getHBondPeakInformation (FILE *msdConfig_file, int nPeaks)
{
	printf("read h bond information from msd.config file...\n");
	fflush (stdout);
	char lineString[1000];
	rewind (msdConfig_file);
	fgets (lineString, 1000, msdConfig_file);

	BOUNDS *peakInfo;
	peakInfo = (BOUNDS *) malloc (nPeaks * sizeof (BOUNDS));

	for (int i = 0; i < nPeaks; ++i)
	{
		fgets (lineString, 1000, msdConfig_file);
		sscanf (lineString, "%f %f\n", &peakInfo[i].lo, &peakInfo[i].hi);
	}

	return peakInfo;
}

float getHBondPeakPosition ()
{
	FILE *hBondThreshold_file;
	hBondThreshold_file = fopen ("hBond.threshold", "r");

	char lineString[1000];
	float thresholdHBondDistance;

	fgets (lineString, 1000, hBondThreshold_file);
	sscanf (lineString, "%f\n", &thresholdHBondDistance);

	return thresholdHBondDistance;
}

void initializeHBondNetworkLogfile ()
{
	system ("rm hBondNetwork_logs/*");
	FILE *hBondNetwork_logfile;
	hBondNetwork_logfile = fopen ("hBondNetwork_logs/hBondNetwork.log", "w");

	// Writing the first line
	fprintf(hBondNetwork_logfile, "# N firstToSecond, N secondToThird, N thirdToFourth, N Mols Zeroth, N Mols First, N Mols Second, N Mols Third, N Mols Fourth\n");

	printf("Created hBond network log file\n");
	fflush (stdout);

	fclose (hBondNetwork_logfile);
}