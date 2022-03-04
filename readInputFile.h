#ifndef READINPUTFILE_H
#define READINPUTFILE_H

DATAFILE_INFO readData (FILE *input, DATA_ATOMS **atoms, DATA_BONDS **bonds, DATA_ANGLES **angles, DATA_DIHEDRALS **dihedrals, DATA_IMPROPERS **impropers);
CONFIG *readConfig (FILE *inputConfigFile, int *nLines_return);
CONFIG *readVWDRadius (FILE *inputVDWConfigFile, int *nLines_return);
DUMPFILE_INFO getDumpFileInfo (FILE *inputDumpFile);

#endif