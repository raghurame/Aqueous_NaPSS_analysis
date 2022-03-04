#ifndef BONDRDF_H
#define BONDRDF_H

void computeBondRDF (DATA_ATOMS *dumpAtoms, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, DATA_BONDS *bonds, CONFIG *inputVectors, DIST_VAR plotVars, int nThreads, float binSize_dist_RDF, float **bondRDF, int *RDFcounter, int currentTimestep);
void printBondRDF (float *bondRDF, int RDFcounter, int nBins_dist_RDF, float binSize_dist_RDF);

#endif