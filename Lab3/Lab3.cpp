#include <iostream>
#include <vector>
#include <omp.h>
#include <ctime>


double analityc_solve(double x, double t, double c) {
	double ksi = x - c * t;
	if (ksi < 0) {
		return 0;
	}
	if (ksi <= 2) {
		return ksi * (2 - ksi);
	}
	return 0;
}


void solve(std::vector<double>& u, double fict_cell, int start, int stop, double c, double h, double tau) {
	double prev_cell, con_cell;
	for (int i = start; i < stop; i++) {
		con_cell = u[i];
		if (i == start) {
			u[i] -= tau * c / h * (u[i] - fict_cell);
		}
		else {
			u[i] -= tau * c / h * (u[i] - prev_cell);
		}
		prev_cell = con_cell;
	}
}


int main(int argc, char*argv[]) {
	int p, N;
	if (argc != 3) {
		std::cerr << "error";
		exit(1);
	}
	else {
		p = std::stoi(argv[1]);
		N = std::stoi(argv[2]);
	}
	
	double h = 10. / (N - 1);
	double t_end = 4;
	double t_start = 0;
	double c = 1;
	double tau = h / c;
	// init vector and initial conditions
	std::vector<double> u(N);
	for (int j = 0; j < N; j++) {
		if (j * h < 0) {
			u[j] = 0;
		}
		else if (j * h <= 2) {
			u[j] = j * h * (2 - j * h);
		}
		else {
			u[j] = 0;
		}
	}

	/*
	double fict_cell = 0;

	while (t_start < t_end) {
		solve(u, &fict_cell, 0, N, c, h, tau);
		t_start += tau;
	}
	std::cout << "X\tu_chisl\tu_analit\n";
	for (int j = 0; j < N; j++) {
		std::cout << j * h << '\t' << u[j] << '\t' << analityc_solve(j * h, t_end, c) << "\n";
	}
	*/
	omp_set_num_threads(p);
	int rank, size, i;
	size = p;
	
	double first = omp_get_wtime();
	#pragma omp parallel shared(u, size, c, h, tau, N, t_end) private(rank, i)
	{	
		int rank = omp_get_thread_num();
		int start, stop;
		start = N / size * rank;
		stop = start + N / size;
		if (N % size) {
			if (rank >= size - N % size) {
				stop += N % size - (size - 1 - rank);
				start += N % size - 1 - (size - 1 - rank);
			}
		}
		printf("Thread number %d: start: %d, stop: %d\n", rank, start, stop);

		double fict_cell;
		double t_start = 0;
		while (t_start < t_end) {
			if (rank == 0) {
				fict_cell = 0;
			}
			if (rank != 0) {
				fict_cell = u[start - 1];
			}
			#pragma omp barrier
			solve(u, fict_cell, start, stop, c, h, tau);
			t_start += tau;
			#pragma omp barrier
		}
	}

	for (int j = 0; j < N; j++) {
		std::cout << j * h << '\t' << u[j] << '\t' << analityc_solve(j * h, t_end, c) << "\n";
	}
	
	std::cout << "Time of working: " << omp_get_wtime() - first << std::endl;

	return 0;
}