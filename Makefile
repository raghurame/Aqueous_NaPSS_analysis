all:
	export OMP_DYNAMIC=true
	gcc -o bondRDF_OOP_entropy bondRDF_OOP_entropy.c -fopenmp -Wall -lm -O3
	gcc -o computeFreeVolume computeFreeVolume.c -fopenmp -Wall -lm -O3
	gcc -o computeMSD computeMSD.c -lm -fopenmp -Wall -O3
	gcc -o bondRDF_ACF bondRDF_ACF.c -lm -fopenmp -Wall -O3
	./bondRDF_OOP_entropy
	./computeFreeVolume
	./computeMSD
	python bondRDF_AUC.py
	python outputMSD.py
