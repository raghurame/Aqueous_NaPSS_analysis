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
#include "generalUtilities.h"
#include "meanSquareDisplacement.h"

void *computeMSD (DATAFILE_INFO datafile, DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, DUMPFILE_INFO dumpfile, int currentDumpstep, DATA_ATOMS **initCoords, MSD_VARS **msdVars, int nPeaks_msd, CONFIG *inputVectors)
{
	float xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo);
	float distance;
	char *outputFilename;
	outputFilename = (char *) malloc (50 * sizeof (char));

	// Store initial coordinates if the currentDumpstep equals 1
	if (currentDumpstep == 1)
	{
		system ("mkdir MSD_logs");
		system ("mkdir MSD_logs/first");
		system ("mkdir MSD_logs/second");
		system ("mkdir MSD_logs/third");
		system ("rm MSD_logs/*");
		system ("rm MSD_logs/first/*");
		system ("rm MSD_logs/second/*");
		system ("rm MSD_logs/third/*");

		// printf("Enter the lower and upper bounds to calculate MSD\n");
		// for (int i = 0; i < nPeaks_msd; ++i)
		// {
		// 	printf("Lower bound for peak %d: ", i + 1); scanf ("%f", &(*msdVars)[i].lowerBound); printf("\n");
		// 	printf("Upper bound for peak %d: ", i + 1); scanf ("%f", &(*msdVars)[i].upperBound); printf("\n");
		// }

		printf("Number of peaks: %d\n", nPeaks_msd);
		fflush (stdout);

		char lineString[1000];

		FILE *msdConfig_file;
		msdConfig_file = fopen ("msd.config", "r");
		fgets (lineString, 1000, msdConfig_file); // Skipping the first line, which contains the number of lines in the config file

		for (int i = 0; i < nPeaks_msd; ++i)
		{
			fgets (lineString, 1000, msdConfig_file);
			sscanf (lineString, "%f %f\n", &(*msdVars)[i].lowerBound, &(*msdVars)[i].upperBound);
		}
		printf("reading peak information from msd.config file...\n");
		fflush (stdout);

		for (int i = 0; i < dumpfile.nAtoms; ++i)
		{
			(*initCoords)[i].id  = dumpAtoms[i].id; (*initCoords)[i].molType = dumpAtoms[i].molType; (*initCoords)[i].atomType = dumpAtoms[i].atomType; (*initCoords)[i].charge = dumpAtoms[i].charge; (*initCoords)[i].x  = dumpAtoms[i].x; (*initCoords)[i].y = dumpAtoms[i].y; (*initCoords)[i].z  = dumpAtoms[i].z; (*initCoords)[i].ix  = dumpAtoms[i].ix; (*initCoords)[i].iy  = dumpAtoms[i].iy; (*initCoords)[i].iz  = dumpAtoms[i].iz; 
		}
	}

	// Calculate mean square displacement if the currentDumpstep is greater than 1
	else if (currentDumpstep > 1)
	{
		for (int i = 0; i < datafile.nBonds; ++i)
		{
			bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType; bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType; bonds[i].x1 = translatePeriodic (dumpAtoms[bonds[i].atom1 - 1].x, dumpAtoms[bonds[i].atom1 - 1].ix, xDist); bonds[i].y1 = translatePeriodic (dumpAtoms[bonds[i].atom1 - 1].y, dumpAtoms[bonds[i].atom1 - 1].iy, yDist); bonds[i].z1 = translatePeriodic (dumpAtoms[bonds[i].atom1 - 1].z, dumpAtoms[bonds[i].atom1 - 1].iz, zDist); bonds[i].x2 = translatePeriodic (dumpAtoms[bonds[i].atom2 - 1].x, dumpAtoms[bonds[i].atom2 - 1].ix, xDist); bonds[i].y2 = translatePeriodic (dumpAtoms[bonds[i].atom2 - 1].y, dumpAtoms[bonds[i].atom2 - 1].iy, yDist); bonds[i].z2 = translatePeriodic (dumpAtoms[bonds[i].atom2 - 1].z, dumpAtoms[bonds[i].atom2 - 1].iz, zDist); bonds[i].xc = (bonds[i].x1 + bonds[i].x2) / 2; bonds[i].yc = (bonds[i].y1 + bonds[i].y2) / 2; bonds[i].zc = (bonds[i].z1 + bonds[i].z2) / 2;

			if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom2Type == inputVectors[0].atom1 && bonds[i].atom1Type == inputVectors[0].atom2))
			{
				for (int j = 0; j < datafile.nBonds; ++j)
				{
					bonds[j].atom1Type = (int) dumpAtoms[bonds[j].atom1 - 1].atomType; bonds[j].atom2Type = (int) dumpAtoms[bonds[j].atom2 - 1].atomType; bonds[j].x1 = translatePeriodic (dumpAtoms[bonds[j].atom1 - 1].x, dumpAtoms[bonds[j].atom1 - 1].ix, xDist); bonds[j].y1 = translatePeriodic (dumpAtoms[bonds[j].atom1 - 1].y, dumpAtoms[bonds[j].atom1 - 1].iy, yDist); bonds[j].z1 = translatePeriodic (dumpAtoms[bonds[j].atom1 - 1].z, dumpAtoms[bonds[j].atom1 - 1].iz, zDist); bonds[j].x2 = translatePeriodic (dumpAtoms[bonds[j].atom2 - 1].x, dumpAtoms[bonds[j].atom2 - 1].ix, xDist); bonds[j].y2 = translatePeriodic (dumpAtoms[bonds[j].atom2 - 1].y, dumpAtoms[bonds[j].atom2 - 1].iy, yDist); bonds[j].z2 = translatePeriodic (dumpAtoms[bonds[j].atom2 - 1].z, dumpAtoms[bonds[j].atom2 - 1].iz, zDist); bonds[j].xc = (bonds[j].x1 + bonds[j].x2) / 2; bonds[j].yc = (bonds[j].y1 + bonds[j].y2) / 2; bonds[j].zc = (bonds[j].z1 + bonds[j].z2) / 2;

					if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom2Type == inputVectors[1].atom1 && bonds[j].atom1Type == inputVectors[1].atom2))
					{
						distance = sqrt (pow ((bonds[j].xc - bonds[i].xc), 2) + pow ((bonds[j].yc - bonds[i].yc), 2) + pow ((bonds[j].zc - bonds[i].zc), 2));

						for (int k = 0; k < nPeaks_msd; ++k)
						{
							if (distance <= (*msdVars)[k].upperBound && distance > (*msdVars)[k].lowerBound)
							{
								FILE *msdOutput;
								snprintf (outputFilename, 50, "MSD_logs/%d-%d_%d-%d_%d.msd", bonds[i].atom1, bonds[i].atom2, bonds[j].atom1, bonds[j].atom2, k + 1);
								msdOutput = fopen (outputFilename, "a");

								fprintf(msdOutput, "%d %f %f %f %f\n", currentDumpstep, distance, bonds[j].xc, bonds[j].yc, bonds[j].zc);

								fclose (msdOutput);
							}
						}
					}
				}
			}
		}
	}

	free (outputFilename);
}