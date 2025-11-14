#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main() {
	int N = 10000;
	std::vector<std::vector<double>> u(N, std::vector<double>(N));

	double manytime = MPI_Wtime();
	# pragma omp parallel for collapse(2)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			u[i][j] = i + j;
		}
	}

	std::cout << "Time: " << MPI_Wtime() - manytime << std::endl;
}