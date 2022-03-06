#ifndef FREEVOLUME_H
#define FREEVOLUME_H

FREEVOLUME_VARS getFreeVolumeVars (DUMPFILE_INFO dumpfile);
int *computeFreeVolume_checkOccupation (int *isOccupied, int i, DUMPFILE_INFO dumpfile, DATA_ATOMS *dumpAtoms, FREEVOLUME_VARS freeVolumeVars, CONFIG *vwdSize, int nThreads);
void computeFreeVolume_getDistribution (int i, int j, FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars, int *isOccupied, DUMPFILE_INFO dumpfile, DATA_ATOMS *dumpAtoms, CONFIG *freeVolumeconfig, int nThreads);
void initializeFreeVolumeDistribution (FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars);
void initializeNBins (FREEVOLUME_VARS *freeVolumeVars);
void resetFreeVolumeDistCounts (FREEVOLUME_DISTRIBUTION **freeVolumeDist, FREEVOLUME_VARS freeVolumeVars);
void printFreeVolumeDistribution (FREEVOLUME_DISTRIBUTION *freeVolumeDist, int atomType, float probeSize, int nBins);
void computeFreeVolume (FREEVOLUME_VARS freeVolumeVars, DATA_ATOMS *dumpAtoms, DUMPFILE_INFO dumpfile, CONFIG *freeVolumeconfig, CONFIG *vwdSize, NLINES_CONFIG entries, int currentDumpstep, int nThreads);

#endif