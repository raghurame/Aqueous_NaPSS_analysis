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
#include "aqueousnapss.h"
#include "fileHandling.h"
#include "readInputFile.h"
#include "generalUtilities.h"
#include "waterOrientation.h"
#include "hBondCorrelation.h"
#include "freeVolume.h"
#include "bondRDF.h"
#include "meanSquareDisplacement.h"

void processLAMMPSTraj (FILE *inputDumpFile, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, DATA_BONDS *bonds, CONFIG *inputVectors, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int nThreads)
{
	// DUMPFILE_INFO dumpfile;
	// dumpfile = getDumpFileInfo (inputDumpFile);

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
	int isTimestep = 0, currentTimestep, currentLine = 0, currentDumpstep = 0, restartDumpStep = 0;

	if (isFileExists ("aqNaPSS.restart"))
	{
		FILE *restartFile;
		restartFile = fopen ("aqNaPSS.restart", "r");
		fgets (lineString, 1000, restartFile);
		sscanf (lineString, "%d", &restartDumpStep);
		printf("\n\nIMPORTANT:\n~~~~~~~~~~\n\nReading restart file. The analysis will begin from %d timeframe in input dump file. If you want to start the analysis from the beginning, then delete the aqNaPSS.restart file and run the program again!\n\n\n");
		fclose (restartFile);
	}

	long int nElements = 0;
	int isNElementsSet = 0;

	// Datafile struct is used to store dump atom information
	DATA_ATOMS *dumpAtoms, *initCoords;
	dumpAtoms = (DATA_ATOMS *) malloc (dumpfile.nAtoms * sizeof (DATA_ATOMS));
	initCoords = (DATA_ATOMS *) malloc (dumpfile.nAtoms * sizeof (DATA_ATOMS));

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

	// MSD variables
	MSD_VARS *msdVars;

	// Checking faults in input dump file
	int currentAtomID = 0, previousAtomID = 0, fault = 0;

	// Limiting the number of timeframes to compute
	int nFreeVolumeCounts = 0, nBondRDFCounts = 0, nWaterOrientationCounts = 0, nHBondComputeCounts = 0, nMSDComputeCounts = 0;

	// Variable to store the number of bonds (dipoles) present in the distance bin
	int *nBonds_inBin_dist_RDF_summation, *nBonds_inBin_dist_RDF_temp;
	nBonds_inBin_dist_RDF_summation = (int *) calloc (nBins_dist_RDF, sizeof (int));

	// The following variables can be used if the input dump file contains scaled coordinates
	double x_temp, y_temp, z_temp;

	// Reading and processing dump information
	while (fgets (lineString, 1000, inputDumpFile) != NULL)
	{
		// Executed at the end of first timeframe
		if (currentLine == (9 + dumpfile.nAtoms) && nElements == 0)
		{
			nElements = getNElements (datafile, dumpAtoms, bonds, inputVectors);
			plotVars.nElements = nElements;
			printf("Allocating memory for %ld elements...\n", nElements);
			allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));
			printf("Memory allocated successfully...\n");

			freeVolumeVars = getFreeVolumeVars (dumpfile);

			// Getting necessary information for mean square displacement calculations
			FILE *msdConfig_file;
			msdConfig_file = fopen ("msd.config", "r");
			printf("Opened msd.config file\n");
			fflush (stdout);
			char lineString2[1000];

			fgets (lineString2, 1000, msdConfig_file);
			sscanf (lineString2, "%d", &nPeaks);

			// Gathering peak information for analysing H-bonds
			peakInfo = getHBondPeakInformation (msdConfig_file, nPeaks);
			peakHBondPosition = getHBondPeakPosition ();
			initializeHBondNetworkLogfile ();
			printf("HBond network file initialized successfully...\n");
			fflush (stdout);

			fclose (msdConfig_file);

			msdVars = (MSD_VARS *) malloc (nPeaks * sizeof (MSD_VARS));
			computeMSD (datafile, dumpAtoms, bonds, dumpfile, currentDumpstep, &initCoords, &msdVars, nPeaks, inputVectors);
		}

		// Main processing loop
		if ((currentDumpstep > restartDumpStep) && (nElements > 0) && (currentLine == 2) && (fault == 0))
		{
			sscanf (lineString, "%d", &currentTimestep);
			printf("Scanning timestep: %d; (%d)...               \n", currentTimestep, currentDumpstep);
			fflush (stdout); 

			if (((currentDumpstep % 1) == 0) && (nBondRDFCounts <= 200))
			{
				nBonds_inBin_dist_RDF_temp = computeBondRDF (dumpAtoms, datafile, dumpfile, bonds, inputVectors, plotVars, nThreads, binSize_dist_RDF, &bondRDF, &RDFcounter, currentTimestep);
				nBondRDFCounts++;

				for (int k = 0; k < nBins_dist_RDF; ++k)
				{
					nBonds_inBin_dist_RDF_summation[k] += nBonds_inBin_dist_RDF_temp[k];
				}
			}

			if (((currentDumpstep % 1) == 0) && (nWaterOrientationCounts <= 200))
			{
				allData_array = computeOrderParameter (allData_array, dumpAtoms, dumpfile, datafile, bonds, inputVectors, currentTimestep, nElements);
				// computeDistribution_OOP (allData_array, plotVars, &distribution_OOP, nThreads);
				computeDistribution_theta (allData_array, plotVars, &distribution_degrees, nThreads);
				nWaterOrientationCounts++;
			}

			// if ((currentDumpstep % 1) == 0 && (nFreeVolumeCounts <= 25))
			// {
			// 	computeFreeVolume (freeVolumeVars, dumpAtoms, dumpfile, freeVolumeconfig, vwdSize, entries, currentDumpstep, nThreads);
			// 	nFreeVolumeCounts++;
			// }

			// if (((currentDumpstep % 5) == 0) && (nHBondComputeCounts <= 400))
			// {
			// 	computeHBonding (dumpAtoms, bonds, datafile, dumpfile, peakInfo, nPeaks, inputVectors, entries, peakHBondPosition, currentDumpstep, nThreads);
			// 	nHBondComputeCounts++;
			// }

			// if (((currentDumpstep % 5) == 0) && (nMSDComputeCounts <= 400))
			// {
			// 	computeMSD (datafile, dumpAtoms, bonds, dumpfile, currentDumpstep, &initCoords, &msdVars, nPeaks, inputVectors);
			// 	nMSDComputeCounts++;
			// }

			isTimestep = 0;
		}
		else if ((nElements > 0) && (currentLine == 2) && (fault == 0))
		{
			sscanf (lineString, "%d", &currentTimestep);
			printf("Skipping timestep: %d; (%d)...               \r", currentTimestep, currentDumpstep);
			fflush (stdout); 
		}

		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			isTimestep = 1;
			currentDumpstep++;
			currentLine = 1;
			isNElementsSet = 0;
			currentAtomID = 0;
			previousAtomID = 0;
			fault = 0;

			// Saving the current progress
			FILE *restartFile;
			restartFile = fopen ("aqNaPSS.restart", "w");
			fprintf(restartFile, "%d", currentDumpstep);
			fclose (restartFile);
		}

		if ((currentLine > 9) && (currentLine < (9 + dumpfile.nAtoms)))
		{
			// Use the following line if the input dump file contains unscaled coordinates as input
			sscanf (lineString, "%d %d %f %f %f %*f %*f %*f %d %d %d\n",
				&dumpAtoms[currentLine - 10].id,
				&dumpAtoms[currentLine - 10].atomType,
				&dumpAtoms[currentLine - 10].x,
				&dumpAtoms[currentLine - 10].y,
				&dumpAtoms[currentLine - 10].z,
				&dumpAtoms[currentLine - 10].ix,
				&dumpAtoms[currentLine - 10].iy,
				&dumpAtoms[currentLine - 10].iz);

			// The following lines can be used if the trajectory dump file contains scaled coordinates
			// Uncomment the following lines if unscaled coordinates are given

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// sscanf (lineString, "%d %d %lf %lf %lf %d %d %d\n",
			// 	&dumpAtoms[currentLine - 10].id,
			// 	&dumpAtoms[currentLine - 10].atomType,
			// 	&x_temp,
			// 	&y_temp,
			// 	&z_temp,
			// 	&dumpAtoms[currentLine - 10].ix,
			// 	&dumpAtoms[currentLine - 10].iy,
			// 	&dumpAtoms[currentLine - 10].iz);

			// x_temp *= (double) xDist;
			// y_temp *= (double) yDist;
			// z_temp *= (double) zDist;

			// dumpAtoms[currentLine - 10].x = (float) x_temp + dumpfile.xlo;
			// dumpAtoms[currentLine - 10].y = (float) y_temp + dumpfile.ylo;
			// dumpAtoms[currentLine - 10].z = (float) z_temp + dumpfile.zlo;

			// x_temp = 0; y_temp = 0; z_temp = 0;

			// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			currentAtomID = dumpAtoms[currentLine - 10].id;

			if (currentAtomID != (previousAtomID + 1))
			{
				printf("Fault found in lammpstrj. This timeframe will be ignored !\n");
				printf("lineString ==> %s\n", lineString);
				printf("previousAtomID: %d; currentAtomID: %d\n", previousAtomID, currentAtomID);
				fault++;
			}
			previousAtomID = currentAtomID;
		}
		currentLine++;
	}

	if (nWaterOrientationCounts > 0)
	{
		// printDistribution_OOP (distribution_OOP, plotVars);
		printDistribution_degrees (distribution_degrees, plotVars);
	}
	if (nBondRDFCounts > 0)
	{
		printBondRDF (bondRDF, RDFcounter, nBins_dist_RDF, binSize_dist_RDF, nBonds_inBin_dist_RDF_summation);
	}

	free (distribution_OOP);
	free (distribution_degrees);
	free (dumpAtoms);
	free (initCoords);
	free (bondRDF);
	free (allData_array);
	free (msdVars);
}