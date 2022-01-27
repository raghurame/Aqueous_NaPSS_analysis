all:
	export OMP_DYNAMIC=true
	gcc -o bondRDF_OOP_entropy bondRDF_OOP_entropy.c -fopenmp -lm -Wall -O3
	./bondRDF_OOP_entropy
