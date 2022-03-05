#ifndef FILEHANDLING_H
#define FILEHANDLING_H

int isFile(const char *name);
int isFileExists (char *inputFilename);
int displayFiles(const char *fileExtension);
char *getInputFileName();
char *getInputFileName_direct (const char *inputFileTemplate);

#endif