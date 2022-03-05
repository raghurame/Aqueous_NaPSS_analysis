#ifndef BONDRDF_ACF_H
#define BONDRDF_ACF_H

void computeACFOfBondRDF (FILE *inputDumpFile);
long long int findNlines (char *logFilename);
void storeLogfileInfo (char *logFilename, ALL_DATA_BONDRDF_ACF **fullData, int nLines, int currentTrajCount);
LOGFILES_VARIABLES_BONDRDF_ACF openLogFiles (FILE *mainDumpfile, ALL_DATA_BONDRDF_ACF **fullData);
void printACF (int i, float *originalDistance, float *acf, int currentTrajCount, float lowerLimit, float upperLimit);
void computeACF (LOGFILES_VARIABLES_BONDRDF_ACF fileVars, ALL_DATA_BONDRDF_ACF *fullData);

#endif