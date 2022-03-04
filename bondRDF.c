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
#include "bondRDF.h"
#include "generalUtilities.h"

void computeBondRDF (DATA_ATOMS *dumpAtoms, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, DATA_BONDS *bonds, CONFIG *inputVectors, DIST_VAR plotVars, int nThreads, float binSize_dist_RDF, float **bondRDF, int *RDFcounter, int currentTimestep)
{
	FILE *bondRDF_logfile;
	char *bondRDF_logfilename;
	bondRDF_logfilename = (char *) malloc (50 * sizeof (char));
	snprintf (bondRDF_logfilename, 50, "bondRDF_logs/%d.rdf", currentTimestep);
	bondRDF_logfile = fopen (bondRDF_logfilename, "w");

	(*RDFcounter)++;

	float xDist = (dumpfile.xhi - dumpfile.xlo), yDist = (dumpfile.yhi - dumpfile.ylo), zDist = (dumpfile.zhi - dumpfile.zlo), simVolume = (xDist * yDist * zDist), bondDensity;
	float xDistHalf = (xDist / 2), yDistHalf = (yDist / 2), zDistHalf = (zDist / 2);
	int nBonds_RDF = 0;

	// Computing bondRDF
	float binStart_dist_RDF = 0, binEnd_dist_RDF, distance;
	int nBins_dist_RDF = (int) (plotVars.maxDist / binSize_dist_RDF);
	int *nBonds_inBin_dist_RDF;
	nBonds_inBin_dist_RDF = (int *) calloc (nBins_dist_RDF, sizeof (int));
	float *nBonds_inBin_dist_RDF_float;
	nBonds_inBin_dist_RDF_float = (float *) malloc (nBins_dist_RDF * sizeof (float));

	float x_translated, y_translated, z_translated;

	omp_set_num_threads (nThreads);
	#pragma omp parallel for
	for (int i = 0; i < datafile.nBonds; ++i)
	{
		if ((bonds[i].atom1Type == inputVectors[0].atom1 && bonds[i].atom2Type == inputVectors[0].atom2) || (bonds[i].atom2Type == inputVectors[0].atom1 && bonds[i].atom1Type == inputVectors[0].atom2))
		{
			bonds[i].atom1Type = (int) dumpAtoms[bonds[i].atom1 - 1].atomType; bonds[i].atom2Type = (int) dumpAtoms[bonds[i].atom2 - 1].atomType;
			bonds[i].x1 = dumpAtoms[bonds[i].atom1 - 1].x; bonds[i].y1 = dumpAtoms[bonds[i].atom1 - 1].y; bonds[i].z1 = dumpAtoms[bonds[i].atom1 - 1].z;
			bonds[i].x2 = dumpAtoms[bonds[i].atom2 - 1].x; bonds[i].y2 = dumpAtoms[bonds[i].atom2 - 1].y; bonds[i].z2 = dumpAtoms[bonds[i].atom2 - 1].z;

			bonds[i].xc = findBondCenter (bonds[i].x1, dumpAtoms[bonds[i].atom1 - 1].ix, bonds[i].x2, dumpAtoms[bonds[i].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
			bonds[i].yc = findBondCenter (bonds[i].y1, dumpAtoms[bonds[i].atom1 - 1].iy, bonds[i].y2, dumpAtoms[bonds[i].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
			bonds[i].zc = findBondCenter (bonds[i].z1, dumpAtoms[bonds[i].atom1 - 1].iz, bonds[i].z2, dumpAtoms[bonds[i].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

			// bonds[i].xc = (bonds[i].x1 + bonds[i].x2) / 2; bonds[i].yc = (bonds[i].y1 + bonds[i].y2) / 2; bonds[i].zc = (bonds[i].z1 + bonds[i].z2) / 2;

			for (int j = 0; j < datafile.nBonds; ++j)
			{
				if ((bonds[j].atom1Type == inputVectors[1].atom1 && bonds[j].atom2Type == inputVectors[1].atom2) || (bonds[j].atom2Type == inputVectors[1].atom1 && bonds[j].atom1Type == inputVectors[1].atom2))
				{
					bonds[j].atom1Type = (int) dumpAtoms[bonds[j].atom1 - 1].atomType; bonds[j].atom2Type = (int) dumpAtoms[bonds[j].atom2 - 1].atomType;
					bonds[j].x1 = dumpAtoms[bonds[j].atom1 - 1].x; bonds[j].y1 = dumpAtoms[bonds[j].atom1 - 1].y; bonds[j].z1 = dumpAtoms[bonds[j].atom1 - 1].z;
					bonds[j].x2 = dumpAtoms[bonds[j].atom2 - 1].x; bonds[j].y2 = dumpAtoms[bonds[j].atom2 - 1].y; bonds[j].z2 = dumpAtoms[bonds[j].atom2 - 1].z;

					bonds[j].xc = findBondCenter (bonds[j].x1, dumpAtoms[bonds[j].atom1 - 1].ix, bonds[j].x2, dumpAtoms[bonds[j].atom2 - 1].ix, dumpfile.xlo, dumpfile.xhi);
					bonds[j].yc = findBondCenter (bonds[j].y1, dumpAtoms[bonds[j].atom1 - 1].iy, bonds[j].y2, dumpAtoms[bonds[j].atom2 - 1].iy, dumpfile.ylo, dumpfile.yhi);
					bonds[j].zc = findBondCenter (bonds[j].z1, dumpAtoms[bonds[j].atom1 - 1].iz, bonds[j].z2, dumpAtoms[bonds[j].atom2 - 1].iz, dumpfile.zlo, dumpfile.zhi);

					// bonds[j].xc = (bonds[j].x1 + bonds[j].x2) / 2; bonds[j].yc = (bonds[j].y1 + bonds[j].y2) / 2; bonds[j].zc = (bonds[j].z1 + bonds[j].z2) / 2;

					x_translated = translatePeriodicDistance (bonds[i].xc, bonds[j].xc, xDistHalf);
					y_translated = translatePeriodicDistance (bonds[i].yc, bonds[j].yc, yDistHalf);
					z_translated = translatePeriodicDistance (bonds[i].zc, bonds[j].zc, zDistHalf);

					// if ((x_translated != bonds[j].xc) || (y_translated != bonds[j].yc) || (z_translated != bonds[j].zc))
					// {
					// 	printf("before translation => [%.3f, %.3f, %.3f] <=> ref: [%.3f, %.3f, %.3f]\nafter translation => [%.3f, %.3f, %.3f] <=> ref: [%.3f, %.3f, %.3f]\n\n", bonds[j].xc, bonds[j].yc, bonds[j].zc, bonds[i].xc, bonds[i].yc, bonds[i].zc, x_translated, y_translated, z_translated, bonds[i].xc, bonds[i].yc, bonds[i].zc);
					// 	sleep (1);
					// }

					// distance = sqrt (pow ((bonds[i].xc - bonds[j].xc), 2) + pow ((bonds[i].yc - bonds[j].yc), 2) + pow ((bonds[i].zc - bonds[j].zc), 2));
					distance = sqrt (pow ((bonds[i].xc - x_translated), 2) + pow ((bonds[i].yc - y_translated), 2) + pow ((bonds[i].zc - z_translated), 2));
					binStart_dist_RDF = 0.0;

					fprintf(bondRDF_logfile, "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n", bonds[i].xc, bonds[i].yc, bonds[i].zc, bonds[j].xc, bonds[j].yc, bonds[j].zc, x_translated, y_translated, z_translated, distance);

					for (int k = 0; k < nBins_dist_RDF; ++k)
					{
						binEnd_dist_RDF = binStart_dist_RDF + binSize_dist_RDF;

						// If i = j, then distance = 0.0
						// neglect the bond pairs where the distance = 0.0
						if (distance > 0 && distance <= binEnd_dist_RDF && distance > binStart_dist_RDF)
							nBonds_inBin_dist_RDF[k]++;

						binStart_dist_RDF = binEnd_dist_RDF;
					}
				}
			}
		}
	}

	for (int i = 0; i < datafile.nBonds; ++i)
	{
		// Finding the bulk density of bond present in inputVectors[1]
		if ((bonds[i].atom1Type == inputVectors[1].atom1 && bonds[i].atom2Type == inputVectors[1].atom2) || (bonds[i].atom2Type == inputVectors[1].atom1 && bonds[i].atom1Type == inputVectors[1].atom2))
			nBonds_RDF++;		
	}

	bondDensity = (float) nBonds_RDF * 20.0 / simVolume;

	binStart_dist_RDF = 0.0;

	for (int i = 0; i < nBins_dist_RDF; ++i)
	{
		binEnd_dist_RDF = binStart_dist_RDF + binSize_dist_RDF;
		nBonds_inBin_dist_RDF_float[i] = (float) nBonds_inBin_dist_RDF[i] / (4.0 * 3.14 * pow (binStart_dist_RDF, 2) * (binEnd_dist_RDF - binStart_dist_RDF));
		(*bondRDF)[i] += (float) nBonds_inBin_dist_RDF_float[i] / bondDensity;
		binStart_dist_RDF = binEnd_dist_RDF;
	}

	fclose (bondRDF_logfile);
}

void printBondRDF (float *bondRDF, int RDFcounter, int nBins_dist_RDF, float binSize_dist_RDF)
{
	FILE *file_bondRDF;
	file_bondRDF = fopen ("bondRDF.output", "w");

	for (int i = 0; i < nBins_dist_RDF; ++i)
	{
		fprintf(file_bondRDF, "%f %f\n", 
			binSize_dist_RDF * (float) i,
			bondRDF[i]/RDFcounter);
	}

	fclose (file_bondRDF);
}