#include <iostream>
#include <mpi.h>
#include <cmath>
#include <vector>


double f(double x) {
	return 4. / (1 + x * x);
}


double integral(int start, int stop, double step) {
	double sum = 0;
	while(start < stop) {
		sum += f(start * step + step / 2.);
		std::cout << start * step + step / 2. << "\n";
		start++;
	}
	return sum * step;
}


int main(int argc, char* argv[]) {
	MPI_Status Status;
	int rank, size;
	double tstart, tend;
	int N;
	if (argc != 2) {
		std::cerr << "error, no input N or too many inputs";
	}
	else {
		N = std::stoi(argv[1]);
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int start, stop;
	double dx = 1. / N;
	if (rank == 0) {
		double ans = integral(0, N, dx);
		std::cout << "One proc ans: " << ans << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	tstart = MPI_Wtime();
	if (rank == 0) {
		stop = N / size;
		for (int i = 1; i < size; i++) {
			start = stop;
			stop = start + N / size;
			if (i < N % size + 1) {
				stop++;
			}
			if (i == size - 1) {
				stop = N;
			}
			MPI_Send(&start, 1, MPI_INT, i, i, MPI_COMM_WORLD);
			MPI_Send(&stop, 1, MPI_INT, i, i, MPI_COMM_WORLD);
		}
		start = 0;
		stop = N / size;
	}

	if (rank != 0) {
		MPI_Recv(&start, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &Status);
		MPI_Recv(&stop, 1, MPI_INT, 0, rank, MPI_COMM_WORLD, &Status);
	}
	
	std::cout << "Rank " << rank << " " << start << " " << stop << '\n';
	double p = integral(start, stop, dx);

	if (rank != 0) {
		MPI_Send(&p, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	
	if (rank == 0) {
		double tmp;	
		std::cout << "I0 " << p << '\n';
		for (int i = 1; i < size; i++) {
			MPI_Recv(&tmp, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Status);
			p += tmp;
			std::cout << "I" << i << " " << tmp << '\n';
		}
		std::cout << "Num_of_procs: " << size << '\n';
		std::cout << "Time of working: " << MPI_Wtime() - tstart << '\n';
		std::cout << "Answer of p procs: " << p << std::endl;	
	}	
	MPI_Finalize();
	return 0;
}
