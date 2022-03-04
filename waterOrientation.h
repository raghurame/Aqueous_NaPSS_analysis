#ifndef WATERORIENTATION_H
#define WATERORIENTATION_H

ORDERPARAMETER *computeOrderParameter (DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, int currentTimestep, unsigned int nElements);
unsigned int getNElements (DATAFILE_INFO datafile, DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, CONFIG *inputVectors);
void computeDistribution_OOP (ORDERPARAMETER *allData_array, DIST_VAR plotVars, DISTRIBUTION **distribution_OOP, int nThreads);
void computeDistribution_theta (ORDERPARAMETER *allData_array, DIST_VAR plotVars, DISTRIBUTION **distribution_degrees, int nThreads);
void printDistribution_OOP (DISTRIBUTION *distribution_OOP, DIST_VAR plotVars);
void printDistribution_degrees (DISTRIBUTION *distribution_degrees, DIST_VAR plotVars);

#endif