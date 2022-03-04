#ifndef MEANSQUAREDISPLACEMENT_H
#define MEANSQUAREDISPLACEMENT_H

void *computeMSD (DATAFILE_INFO datafile, DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, DUMPFILE_INFO dumpfile, int currentDumpstep, DATA_ATOMS **initCoords, MSD_VARS **msdVars, int nPeaks_msd, CONFIG *inputVectors);

#endif