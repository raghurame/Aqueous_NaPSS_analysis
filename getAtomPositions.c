#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

int main(int argc, char const *argv[])
{
	if (argc < 4)
	{
		printf("Required args:\n~~~~~~~~~~~~~~\n\nargv[0] = ./program\nargv[1] = input dump file\nargv[2] = required atom index to capture from input dump\nargv[3] = output file to dump positions.\n\n");
		exit (1);
	}

	FILE *inputDump, *outputPositions;
	inputDump = fopen (argv[1], "r");
	outputPositions = fopen (argv[3], "w");

	char lineString[1000];
	int requiredAtomIndex = atoi (argv[2]), atomIndex;

	while (fgets (lineString, 1000, inputDump) != NULL)
	{
		sscanf (lineString, "%d\n", &atomIndex);
		if (atomIndex == requiredAtomIndex)
		{
			fprintf(outputPositions, "%s", lineString);
		}
	}

	fclose (inputDump);
	return 0;
}