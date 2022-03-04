all:
	export OMP_DYNAMIC=true
	gcc -c aqueousnapss.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c fileHandling.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c readInputFile.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c generalUtilities.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c waterOrientation.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c hBondCorrelation.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c freeVolume.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c bondRDF.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c meanSquareDisplacement.c -Wall -fstack-protector -g -fopenmp -lm
	gcc -c main.c -Wall -fstack-protector -g -fopenmp -lm
	gcc main.o aqueousnapss.o fileHandling.o readInputFile.o generalUtilities.o waterOrientation.o hBondCorrelation.o freeVolume.o bondRDF.o meanSquareDisplacement.o -o main -lm -fopenmp -Wall -g -fstack-protector
