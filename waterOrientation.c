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
#include "waterOrientation.h"
#include "generalUtilities.h"

ORDERPARAMETER *computeOrderParameter (ORDERPARAMETER *allData_array, DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, int currentTimestep, unsigned int nElements)
{
	// FILE *allData;
	// char *allData_string;
	// allData_string = (char *) malloc (50 * sizeof (char));
	// sprintf (allData_string, "logs/allData_%d.oop", currentTimestep);
	// allData = fopen (allData_string, "w");

	// fprintf(allData, "atom1, atom2, atom3, atom4, distance, angle (rad), angle (deg), OOP\n");

	// ORDERPARAMETER *allData_array;
	// allData_array = (ORDERPARAMETER *) malloc (nElements * sizeof (ORDERPARAMETER));

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
	// free (allData_string);

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
				if ((allData_array[k].theta_deg <= currentBounds.binEnd_deg) && (allData_array[k].theta_deg > currentBounds.binStart_deg) && (allData_array[k].distance <= currentBounds.binEnd_dist) && (allData_array[k].distance > currentBounds.binStart_dist))
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
	fprintf(file_distribution_degrees_info, "binStart_dist: %f\nbinEnd_dist: %f\nbinStart_deg: %f\nbinEnd_deg: %f\nnBins_dist: %d\nnBins_deg: %d\nsize_degrees: %d\n\n",
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