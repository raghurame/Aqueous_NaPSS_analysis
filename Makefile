all:
	export OMP_DYNAMIC=true
	gcc -o bondRDF_OOP_entropy bondRDF_OOP_entropy.c -fopenmp -Wall -lm -O3
	./bondRDF_OOP_entropy
