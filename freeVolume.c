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
#include "freeVolume.h"
#include "generalUtilities.h"

FREEVOLUME_VARS getFreeVolumeVars (DUMPFILE_INFO dumpfile)
{
	FREEVOLUME_VARS freeVolumeVars;

	FILE *freeVolumeVars_file;
	freeVolumeVars_file = fopen ("freeVolumeVars.config", "r");

	char lineString[1000];

	fgets (lineString, 1000, freeVolumeVars_file);
	sscanf (lineString, "%f", &freeVolumeVars.minProbeSize);

	fgets (lineString, 1000, freeVolumeVars_file);
	sscanf (lineString, "%f", &freeVolumeVars.maxProbeSize);

	fgets (lineString, 1000, freeVolumeVars_file);
	sscanf (lineString, "%f", &freeVolumeVars.delProbeSize);

	fgets (lineString, 1000, freeVolumeVars_file);
	sscanf (lineString, "%f", &freeVolumeVars.binEnd_dist);

	fgets (lineString, 1000, freeVolumeVars_file);
	sscanf (lineString, "%f", &freeVolumeVars.binSize_dist);

	// // Variables regarding probe size
	// printf("Enter the minimum size of probe:\t");
	// scanf ("%f", &freeVolumeVars.minProbeSize);
	// printf("Enter the maximum size of probe:\t");
	// scanf ("%f", &freeVolumeVars.maxProbeSize);
	// printf("Enter the del (size) for probe:\t");
	// scanf ("%f", &freeVolumeVars.delProbeSize);

	freeVolumeVars.nBins_probeSweep = (int) ((freeVolumeVars.maxProbeSize - freeVolumeVars.minProbeSize) / freeVolumeVars.delProbeSize) + 1;

	// // Variables regarding 3D probe movement
	// // delDistance can usually be set at 1/4th of probe radius
	freeVolumeVars.xLength = (dumpfile.xhi - dumpfile.xlo); freeVolumeVars.yLength = (dumpfile.yhi - dumpfile.ylo); freeVolumeVars.zLength = (dumpfile.zhi - dumpfile.zlo);

	// // Variables for calculating free volume distribution
	freeVolumeVars.binStart_dist = 0;
	// printf("Free volume distribution will be calculated from the center of respective atoms\n\n");
	// printf("Enter the max. distance to consider for free volume distribution (float type):\t"); scanf ("%f", &freeVolumeVars.binEnd_dist);
	// printf("Enter the bin size (number of bins will be calculated based on max. distance and bin size):\t"); scanf ("%f", &freeVolumeVars.binSize_dist);
	freeVolumeVars.nBins_dist = (int) ((freeVolumeVars.binEnd_dist - freeVolumeVars.binStart_dist) / freeVolumeVars.binSize_dist) + 1;

	printf("Successfully read free volume variables from config file...\n");
	fflush (stdout);

	return freeVolumeVars;
}

int *computeFreeVolume_checkOccupation (int *isOccupied, int i, DUMPFILE_INFO dumpfile, DATA_ATOMS *dumpAtoms, FREEVOLUME_VARS freeVolumeVars, CONFIG *vwdSize, int nThreads)
{
	DATA_ATOMS probePosition;
	float distance, x1, y1, z1, x2, y2, z2;
	// int *isOccupied;
	long long int arraySize = (long long int)(freeVolumeVars.nBins_dist_x + 1) * (long long int)(freeVolumeVars.nBins_dist_y + 1) * (long long int)(freeVolumeVars.nBins_dist_z + 1), index1d, progress = 0;
	float progressPercent;
	// isOccupied = (int *) calloc (arraySize, sizeof (int));

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
						exit (1);
					}

					if (distance < (freeVolumeVars.currentProbeSize + (vwdSize[dumpAtoms[m].atomType - 1].radius / 2)))
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
	long long int index1d, progress = 0;
	// float progressPercent;
	// isOccupied = (int *) malloc (arraySize * sizeof (int));

	float xDistHalf = (dumpfile.xhi - dumpfile.xlo) / 2, yDistHalf = (dumpfile.yhi - dumpfile.ylo) / 2, zDistHalf = (dumpfile.zhi - dumpfile.zlo) / 2;

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

						x2 = translatePeriodicDistance (x1, x2, xDistHalf);
						y2 = translatePeriodicDistance (y1, y2, yDistHalf);
						z2 = translatePeriodicDistance (z1, z2, zDistHalf);

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

void computeFreeVolume (FREEVOLUME_VARS freeVolumeVars, DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int currentDumpstep, int nThreads)
{
	DATA_ATOMS probePosition;
	freeVolumeVars.currentProbeSize = freeVolumeVars.minProbeSize;

	FREEVOLUME_DISTRIBUTION *freeVolumeDist;
	freeVolumeDist = (FREEVOLUME_DISTRIBUTION *) malloc (freeVolumeVars.nBins_dist * sizeof (FREEVOLUME_DISTRIBUTION));

	initializeFreeVolumeDistribution (&freeVolumeDist, freeVolumeVars);

	long long int arraySize;

	omp_set_num_threads (nThreads);

	// Probe size is set. It'll vary in a loop, from minimum to maximum size
	for (int i = 0; i < freeVolumeVars.nBins_probeSweep; ++i)
	{
		initializeNBins (&freeVolumeVars);

		probePosition.x = dumpfile.xlo; probePosition.y = dumpfile.ylo; probePosition.z = dumpfile.zlo;

		int *isOccupied;
		arraySize = (long long int)freeVolumeVars.nBins_dist_x * (long long int)freeVolumeVars.nBins_dist_y * (long long int)freeVolumeVars.nBins_dist_z;
		isOccupied = (int *) malloc (arraySize * sizeof (int));
		isOccupied = computeFreeVolume_checkOccupation (isOccupied, i, dumpfile, dumpAtoms, freeVolumeVars, vwdSize, nThreads);

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