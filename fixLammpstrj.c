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
#include "structDefinitions.h"
#include "fixLammpstrj.h"
#include "readInputFile.h"

void fixLammpstrj (const char *inputDumpFilename)
{
	FILE *inputDump;
	inputDump = fopen (inputDumpFilename, "r");

	DUMPFILE_INFO dumpfile;
	dumpfile = getDumpFileInfo (inputDump);

	char lineString[1000], **fullDumpText;
	fullDumpText = (char **) malloc ((dumpfile.nAtoms + 9) * sizeof (char *));

	int nTimeframes, currentAtomID = 0, previousAtomID = 0, fault = 0;
	nTimeframes = countNTimeframes (inputDump);

	printf("Checking the input lammpstrj file for write errors...\n");
	fflush (stdout);

	for (int i = 0; i < nTimeframes; ++i)
	{
		fullDumpText[i] = (char *) malloc (1500 * sizeof (char));
	}

	rewind (inputDump);

	while ((fgets (lineString, 1000, inputDump) != NULL))
	{
		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			currentAtomID = 0;
			previousAtomID = 0;
			fault = 0;
			snprintf (fullDumpText[0], 1500, "ITEM: TIMESTEP");

			for (int j = 0; j < 8; ++j)
			{
				fgets (lineString, 1000, inputDump);
				snprintf (fullDumpText[j + 1], 1500, "%s", lineString);
			}

			for (int j = 0; j < dumpfile.nAtoms; ++j)
			{
				fgets (lineString, 1000, inputDump);
				sscanf (lineString, "%d\n", &currentAtomID);
				snprintf (fullDumpText[j + 9], 1500, "%s", lineString);

				if (currentAtomID != (previousAtomID + 1))
					fault++;

				previousAtomID = currentAtomID;
			}

			if (fault == 0)
			{
				for (int i = 0; i < (dumpfile.nAtoms + 9); ++i)
				{
					
				}
			}
		}
	}

	fclose (inputDump);
}