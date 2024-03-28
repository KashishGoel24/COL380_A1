omp: omp_pointers.cpp
	g++-13 omp_pointers.cpp -fopenmp -o ./omp

pthread: pth_pointers.cpp
	g++-13 pth_pointers.cpp -o ./pth -lpthread

serial: serial.cpp
	g++-13 serial.cpp -o ./seq

clean:
	rm -f ./pth ./seq ./omp