#ifndef AQUEOUSNAPSS_H
#define AQUEOUSNAPSS_H

void processLAMMPSTraj (FILE *inputDumpFile, DATAFILE_INFO datafile, DATA_BONDS *bonds, CONFIG *inputVectors, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int nThreads);

#endif