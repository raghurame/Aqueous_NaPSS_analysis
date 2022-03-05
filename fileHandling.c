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

int isFile(const char *name)
{
	DIR *directory = opendir (name);
	if (directory!=NULL)
	{
		closedir(directory);
		return 0;
	}
	if(errno==ENOTDIR)
	{
		return 1;
	}

	return -1;
}

int isFileExists (char *inputFilename)
{
	FILE *checking;

	if (checking = fopen (inputFilename, "r"))
	{
		fclose (checking);
		return 1;
	}
	return 0;
}

int displayFiles(const char *fileExtension)
{
	int nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;
	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			printf("%d --> %s\n", nFiles, filePointer -> d_name);
		}
	}
	return nFiles;
}

char *getInputFileName()
{
	int nFiles = 0;
	char *inputFileName, fileExtension[200];
	int fileRequired;
	inputFileName = (char *) malloc (200 * sizeof (char));

	getFilenameAgain:
	printf("Enter the file extension or a match string to search in current directory...\n --> "); scanf ("%s", fileExtension);
	printf("\n");
	nFiles = displayFiles (fileExtension);

	if (nFiles > 0)
	{
		fprintf(stdout, "\nWhich file would you like to input? Enter a number between (1 to %d): ", nFiles); 
		fflush (stdout);
		scanf ("%d", &fileRequired);
	}
	else
	{
		printf("No files found with the match string. Try again!\n"); goto getFilenameAgain;
	}

	nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;

	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			if (fileRequired == nFiles)
			{
				strcpy (inputFileName, filePointer -> d_name);
			}
		}
	}
	return inputFileName;
}

char *getInputFileName_direct (const char *inputFileTemplate)
{
	int nFiles = 0;
	char *inputFileName;
	inputFileName = (char *) malloc (100 * sizeof (char));

	nFiles = displayFiles (inputFileTemplate);

	if (nFiles == 1)
	{
		DIR *parentDirectory;
		parentDirectory = opendir ("./");

		struct dirent *filePointer;

		/* Scan all the files using filePointer */
		while ((filePointer = readdir (parentDirectory)))
		{
			if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, inputFileTemplate))
			{
				strcpy (inputFileName, filePointer -> d_name);
			}
		}

		return inputFileName;
	}
	else if (nFiles > 1)
	{
		printf("\nERROR:\n~~~~~~\n\nMore than one file found with *.%s extension. Args must be passed !\n\nREQUIRED ARGS:\n~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = dump filename\nargv[2] = data filename\n\n", inputFileTemplate);
		exit (1);
		snprintf (inputFileName, 100, "NULL");
		return inputFileName;
	}
	else if (nFiles == 0)
	{
		printf("\nERROR:\n~~~~~~\n\nNo file found with *.%s extension. Args must be passed !\n\nREQUIRED ARGS:\n~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = dump filename\nargv[2] = data filename\n\n", inputFileTemplate);
		exit (1);
	}
}