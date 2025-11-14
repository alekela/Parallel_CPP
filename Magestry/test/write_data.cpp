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


namespace fs = std::filesystem;


void write_data_parallel(std::vector<std::vector<double>> u, std::vector<std::vector<double>> v, double time, double Re, double dx, double dy,
					int start_x, int stop_x, int start_y, int stop_y, 
					int fict, int rank_x, int rank_y, std::string filename, int size) {
    if (rank_x + rank_y == 0) {
        if (!fs::exists(filename)) {
            fs::create_directory(filename);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string outpath = "Test_parallel.csv";
    std :: ostringstream buffer;
    
    if (rank_x + rank_y == 0) {
        buffer << "X,Y,U_exp,V_exp,u_theory,v_theory\n";
    }
	double x, y;
    for (int i = fict; i < stop_y - start_y + fict; i++) {
        for (int j = fict; j < stop_x - start_x + fict; j++) {
			x = (start_x + (j - fict)) * dx + dx / 2.;
			y = (start_y + (i - fict)) * dy + dy / 2.;
        	buffer << x << "," << y << "," << u[i][j] << "," << v[i][j] << "," << "\n";
        }
    }
    std::string local_data = buffer.str();
    int local_size = local_data.size();

    int offset = 0;
    MPI_Exscan(&local_size, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, outpath.c_str(),MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_at_all(fh, offset, local_data.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}


void write_data(std::vector<std::vector<double>> u, std::vector<std::vector<double>> v, double time, double Re, double dx, double dy,
					int start_x, int stop_x, int start_y, int stop_y, 
					int fict, int rank_x, int rank_y, std::string filename, int size, MPI_Status* Status) {
    if (rank_x + rank_y == 0) {
        if (!fs::exists(filename)) {
            fs::create_directory(filename);
        }
    }
	std::string outpath = "Test_normal.csv";
    	std :: ostringstream buffer;
    
    	if (rank_x + rank_y == 0) {
     		buffer << "X,Y,U_exp,V_exp,u_theory,v_theory\n";
    	}

	double x, y;
    	for (int i = fict; i < stop_y - start_y + fict; i++) {
        	for (int j = fict; j < stop_x - start_x + fict; j++) {
			x = (start_x + (j - fict)) * dx + dx / 2.;
			y = (start_y + (i - fict)) * dy + dy / 2.;
        		buffer << x << "," << y << "," << u[i][j] << "," << v[i][j] << "," << "\n";
        	}
    	}
    	std::string local_data = buffer.str();
    	int local_size = local_data.size();
	
	if (rank_x + rank_y != 0){
		MPI_Send(&local_size, 1, MPI_INT, 0, size, MPI_COMM_WORLD);
	}
	std::vector<int> lengths(size);
	if (rank_x + rank_y == 0) {
		for (int i = 1; i < size; i++) {
			MPI_Recv(&(lengths[i]), 1, MPI_INT, i, size, MPI_COMM_WORLD, Status);
		}
	}
	if (rank_x + rank_y != 0){
		MPI_Send(&(local_data[0]), local_size, MPI_CHAR, 0, size, MPI_COMM_WORLD);
	}

	if (rank_x + rank_y == 0) {
		std::ofstream out;
    		out.open(outpath); 
		for (int i = 1; i < size; i++) {
			std::vector<char> buffwrite(lengths[i]);
			MPI_Recv(buffwrite.data(), lengths[i], MPI_CHAR, i, size, MPI_COMM_WORLD, Status);
			for (const auto &e : buffwrite) out << e;
		}
		out.close();
	}
}


int main(int argc, char* argv[]) {
	int rank, size;
	int N;
	if (argc != 2) {
		std::cerr << "error, no input N or too many inputs";
	}
	else {
		N = std::stoi(argv[1]);
	}
	int ddt = N / 100;
	
	double t_start = 0;
	double t_end = 0.03 / ddt / ddt;
	int Nx = N;
	int Ny = N;
	double dt = 0.001 / ddt / ddt;
	double Re = 100;
	int fict = 1;
	double dx = 1. / Nx;
	double dy = 1. / Ny;
	std::string message;
	std::string filename = "Out_csvs_N_" + std::to_string(N);
	int iterwrite = 10000;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status Status;
	MPI_Request Request;

	filename += "_p_" + std::to_string(size);

	int px, py;
	for (int i = sqrt(size); i > 0; i--) {
		if (size % i == 0) {
			py = i;
			px = size / i;
			break;
		}
	}
	// std::cout << px << " " << py << std::endl;
	int start_x, stop_x, start_y, stop_y;
	int rank_x, rank_y;
	rank_x = rank % px;
	rank_y = rank / px;

	start_x = Nx / px * rank_x;
	stop_x = start_x + Nx / px;
	if (Nx % px) {
		if (rank_x >= px - Nx % px) {
			stop_x += Nx % px - (px - 1 - rank_x);
			start_x += Nx % px - 1 - (px - 1 - rank_x);
		}
	}

	start_y = Ny / py * rank_y;
	stop_y = start_y + Ny / py;
	if (Ny % py) {
		if (rank_y >= py - Ny % py) {
			stop_y += Ny % py - (py - 1 - rank_y);
			start_y += Ny % py - 1 - (py - 1 - rank_y);
		}
	}

	int neighbors[4];
	neighbors[2] = rank_y == 0 ? -1 : px * (rank_y - 1) + rank_x;
	neighbors[1] = rank_x == 0 ? -1 : px * rank_y + rank_x - 1;
	neighbors[0] = rank_y == py - 1 ? -1 : px * (rank_y + 1) + rank_x;
	neighbors[3] = rank_x == px - 1 ? -1 : px * rank_y + rank_x + 1;

	// message = "Rank " + std::to_string(rank) + " start_x: " + std::to_string(start_x) + " stop_x: " + std::to_string(stop_x) + " start_y: " + std::to_string(start_y) + + " stop_y: " + std::to_string(stop_y);
	// std::cout << message << std::endl;

	std::vector<std::vector<double>> u(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	std::vector<std::vector<double>> v(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));

	double manytime = MPI_Wtime();
	write_data_parallel(u, v, t_start, Re, dx, dy, start_x, stop_x, start_y, stop_y, fict, rank_x, rank_y, filename, size);
	if (rank == 0) {
		std::cout << "\nTime of parallel working: " << MPI_Wtime() - manytime << "\n";
	}
	double onetime = MPI_Wtime();
	write_data(u, v, t_start, Re, dx, dy, start_x, stop_x, start_y, stop_y, fict, rank_x, rank_y, filename, size, &Status);
	if (rank == 0) {
		std::cout << "\nTime of normal working: " << MPI_Wtime() - onetime << "\n";
	}
	
	MPI_Finalize();
	return 0;
}



