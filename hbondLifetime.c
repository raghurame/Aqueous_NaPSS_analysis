#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

/*
 * ARGS TO PASS:
 * ~~~~~~~~~~~~
 *
 * argv[0] = ./program
 * argv[1] = input main traj file
 *
 */

int getNextTimestep (FILE *inputTraj)
{
	int currentTimestep, currentLine = 0;
	char lineString[1000];

	while ((fgets (lineString, 1000, inputTraj) != NULL))
	{
		if (strstr (lineString, "ITEM: TIMESTEP"))
		{
			printf("%s\n", lineString);
			fgets (lineString, 1000, inputTraj);
			sscanf (lineString, "%d", &currentTimestep);
			return currentTimestep;
		}
	}

	return 0;
}

int main(int argc, char const *argv[])
{
	FILE *inputTraj;
	inputTraj = fopen (argv[1], "r");

	int currentTimestep = 1;

	do
	{
		currentTimestep = getNextTimestep (inputTraj);
		printf("currentTimestep: %d\n", currentTimestep);
		usleep (1000000);
	} while (currentTimestep > 0);

	return 0;
}