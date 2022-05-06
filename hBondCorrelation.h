#ifndef HBONDCORRELATION_H
#define HBONDCORRELATION_H

DATA_ATOMS *assignPeaks (DATA_ATOMS *dumpAtoms, DATA_ATOMS *dumpAtomsMod, DATA_BONDS *bonds, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, BOUNDS *peakInfo, int nPeaks, CONFIG *inputVectors, int nThreads);
void analyzeHBondNetwork (DATA_ATOMS *dumpAtomsMod, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, CONFIG *inputVectors, BOUNDS *peakInfo, int nPeaks, float peakHBondPosition, int currentDumpstep, int nThreads);
float *findHBondCorrelation (const char *filename, int nTimeframes, int nThreads);
void calculateSumCorrelation (float *correlation, float **sumCorrelation, int nTimeframes);
float *calculateAverageCorrelation (float *sumCorrelation, int nFiles, int nTimeframes);
void computeHBondCorrelation2 (const char *fileTemplate, int nTimeframes, int nThreads, int nFiles_0);
void computeHBondCorrelation (FILE *inputDumpFile, int nThreads);
void computeHBonding (DATA_ATOMS *dumpAtoms, DATA_BONDS *bonds, DATAFILE_INFO datafile, DUMPFILE_INFO dumpfile, BOUNDS *peakInfo, int nPeaks, CONFIG *inputVectors, NLINES_CONFIG entries, float peakHBondPosition, int currentDumpstep, int nThreads);
BOUNDS *getHBondPeakInformation (FILE *msdConfig_file, int nPeaks);
float getHBondPeakPosition ();
void initializeHBondNetworkLogfile ();
int countNFiles (const char *fileTemplate);

#endif