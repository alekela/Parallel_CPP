#include <iostream>
#include <mpi.h>
#include <vector>
#include <cmath>
#include <fstream>


double analityc_solve(double k, double L, double u0, double t, double x, double prec) {
	int m = (int) std::sqrt(-std::log(prec) * L * L / k / M_PI / M_PI / t) / 2 + 1;
	double sum = 0;
	for (int i = 0; i <= m; i++) {
		sum += std::exp(-k * M_PI * M_PI * (2 * i + 1) * (2 * i + 1) * t / L / L) / (2 * i + 1) * std::sin(M_PI * (2 * i + 1) * x / L);
	}
	return sum * 4 * u0 / M_PI;
}


void solve(std::vector<double>* T, int size, double k, double h, double tau) {
	std::vector<double> new_T(size - 2);
	for (int i = 1; i < size - 1; i++) {
		new_T[i - 1] = tau * k / h / h * ((*T)[i + 1] - 2 * (*T)[i] + (*T)[i - 1]) + (*T)[i];
	}
	for (int i = 1; i < size - 1; i++) {
		(*T)[i] = new_T[i - 1];
	}
}


void sending_p(std::vector<double>* T, int n, int rank, int size, MPI_Status* Status) {
	if (rank != size - 1) {
		MPI_Send(&((*T)[n]), 1, MPI_DOUBLE, rank +1, rank +1, MPI_COMM_WORLD);
	}
	if (rank != 0) {
		MPI_Recv(&((*T)[0]), 1, MPI_DOUBLE, rank -1, rank , MPI_COMM_WORLD, Status);
	}
	if (rank != 0) {
		MPI_Send(&((*T)[1]), 1, MPI_DOUBLE, rank -1, rank -1, MPI_COMM_WORLD);
	}
	if (rank != size - 1) {
		MPI_Recv(&((*T)[n + 1]), 1, MPI_DOUBLE, rank +1, rank , MPI_COMM_WORLD, Status);
	}
}


void sending_1(std::vector<double>* T, int n, int rank, int size, MPI_Status* Status) {
	if (rank % 2 != 0) {
		if (rank != 0) {
			MPI_Send(&((*T)[1]), 1, MPI_DOUBLE, rank -1, rank -1, MPI_COMM_WORLD);
			MPI_Recv(&((*T)[0]), 1, MPI_DOUBLE, rank -1, rank , MPI_COMM_WORLD, Status);
		}
		if (rank != size - 1) {
			MPI_Recv(&((*T)[n + 1]), 1, MPI_DOUBLE, rank +1, rank , MPI_COMM_WORLD, Status);
			MPI_Send(&((*T)[n]), 1, MPI_DOUBLE, rank +1, rank +1, MPI_COMM_WORLD);
		}
	}
	if (rank % 2 == 0) {
		if (rank != size - 1) {
			MPI_Recv(&((*T)[n + 1]), 1, MPI_DOUBLE, rank +1, rank , MPI_COMM_WORLD, Status);
			MPI_Send(&((*T)[n]), 1, MPI_DOUBLE, rank +1, rank +1, MPI_COMM_WORLD);
		}
		if (rank != 0) {
			MPI_Send(&((*T)[1]), 1, MPI_DOUBLE, rank -1, rank -1, MPI_COMM_WORLD);
			MPI_Recv(&((*T)[0]), 1, MPI_DOUBLE, rank -1, rank , MPI_COMM_WORLD, Status);
		}
	}
}


int main(int argc, char* argv[]) {
	double xstart = 0;
	double xend = 1;
	int N;
	if (argc != 2) {
		std::cerr << "error, no input N or too many inputs";
	}
	else {
		N = std::stoi(argv[1]);
	}
	double k = 1;
	double tend = 0.0001;
	double h = (xend - xstart) / N;

	double tau = 0.5 * h * h / k;
	int Nt = (int) (tend / tau); 

	double u0 = 1;
	double uleft = 0;
	double uright = 0;
	double points[11] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7001, 0.8, 0.9, 1.};

	MPI_Status Status;
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double onetime = MPI_Wtime();
		std::vector<double> T0(N + 2);
		std::vector<double> T_theory(N);
		for (int i = 1; i < N + 1; i++) {
			T0[i] = u0;
		}
		T0[0] = uleft;
		T0[N + 1] = uright;
		for (int i = 0; i < Nt; i++) {
			solve(&T0, N + 2, k, h, tau);
		}
		std::cout << "One proc time: " << MPI_Wtime() - onetime << "\n";
		// saving analitic solution
		/*
    	std::ofstream outfile("T_analit(X).csv");
    	std::string tmp;
   		tmp = "Time,X,T_chisl,T_analit\n";
   		outfile << tmp;
		int index = 0;
    	for (int i = 0; i < 11; i++) {
       		tmp = "";
        	tmp += std::to_string(tau*Nt) + ',';
       		tmp += std::to_string(xstart + 0.1 * (i)) + ',';
			index = points[i] / h;
			if (i >= 5) {
				index++;
			}
        	tmp += std::to_string(T0[index]) + ',';
        	tmp += std::to_string(analityc_solve(k, xend - xstart, u0, tau*Nt, xstart + 0.1 * (i), 0.00001)) + '\n';
       		outfile << tmp;
    	}
    	outfile.close();
		*/
		int index = 0;
		std::cout << "T_chisl: ";
    	for (int i = 0; i < 11; i++) {
			index = points[i] / h;
			if (i >= 5) {
				index++;
			}
        	std::cout << T0[index] << ',';
		}
		std::cout << "\nT_analit: ";
		for (int i = 0; i < 11; i++) {
        	std::cout << analityc_solve(k, xend - xstart, u0, tau*Nt, xstart + 0.1 * (i), 0.00001) << ',';
    	}
		std::cout << '\n';
	}

	MPI_Barrier(MPI_COMM_WORLD);
	// parallel part
	double manytime = MPI_Wtime();
	int start, stop;
	int conds[2];
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
			conds[0] = start;
			conds[1] = stop;
			MPI_Send(conds, 2, MPI_INT, i, i, MPI_COMM_WORLD);
		}
		start = 0;
		stop = N / size;
	}
	
	if (rank != 0) {
		MPI_Recv(conds, 2, MPI_INT, 0, rank, MPI_COMM_WORLD, &Status);
		start = conds[0];
		stop  = conds[1];
	}
	
	// sync before starting to solve
	MPI_Barrier(MPI_COMM_WORLD);
	std::vector<double> T(stop - start + 2);
	for (int i = 1; i < stop - start + 1; i++) {
		T[i] = 1;
	}
	for (int i = 0; i < Nt; i++) {
		if (rank == 0) {
			T[0] = uleft;
		}
		if (rank == size - 1) {
			T[stop - start + 1] = uright;
		}
		sending_p(&T, stop - start, rank, size, &Status);

		solve(&T, stop - start + 2, k, h, tau);
	}
	int n = stop - start;
	/*
	if (rank != 0) {
		MPI_Send(&n, 1, MPI_INT, 0, rank + size, MPI_COMM_WORLD);
		MPI_Send(&(T[1]), n + 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}*/
	if (rank == 0) {
		std::cout << "Parallel results:" << std::endl;
		/*int counter = 0;
		int ind_count = 0;
		int indexes[11];
    	for (int i = 0; i < 11; i++) {
			indexes[i] = i * N / 10;
			if (i >= 5) {
				indexes[i]++;
			}
		}
		
		std::cout << "T_chisl: ";
		for (int j = 0; j < n + 1 + (size == 1); j++) {
			if (T[j] != T[j]) {
				std::cout << "Nan detected!" << std::endl;
				break;
			}

			if (counter == indexes[ind_count]) {
				std::cout << T[j] << ',';
				ind_count++;
			}
			counter++;
		}
		for (int i = 1; i < size; i++) {
			MPI_Recv(&n, 1, MPI_INT, i, i + size, MPI_COMM_WORLD, &Status);
			MPI_Recv(&(T[1]), n + 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &Status);
			for (int j = 1; j < n + 1 + (i == size - 1); j++) {
				if (T[j] != T[j]) {
					std::cout << "Nan detected!" << std::endl;
					break;
				}

				if (counter == indexes[ind_count]) {
					std::cout << T[j] << ',';
					ind_count++;
				}
				counter++;
			}
		}*/
		std::cout << "\nTime of working: " << MPI_Wtime() - manytime << "\n";
	}
	MPI_Finalize();
		
}