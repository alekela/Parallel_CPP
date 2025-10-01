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


double u_theory(double x, double y, double t, double Re) {
	return 0.75 - 0.25 * (1. / (1 + exp((-4. * x + 4. * y - t) * Re / 32.)));
}


double v_theory(double x, double y, double t, double Re) {
	return 0.75 + 0.25 * (1. / (1 + exp((-4. * x + 4. * y - t) * Re / 32.)));
}


void InitialState(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, 
					int start_x, int stop_x, int start_y, int stop_y, 
					int fict, double dx, double dy, double Re) {
	for (int i = start_y - fict; i < stop_y + fict; i++) {
		for (int j = start_x - fict; j < stop_x + fict; j++) {
			u[i - start_y + fict][j - start_x + fict] = u_theory(j * dx + dx / 2., i * dy + dy / 2., 0, Re);
			v[i - start_y + fict][j - start_x + fict] = v_theory(j * dx + dx / 2., i * dy + dy / 2., 0, Re);
		}
	}
}


void Boundary(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& v, double time, double Re, double dx, double dy,
					int start_x, int stop_x, int start_y, int stop_y, 
					int fict, int rank_x, int rank_y, int px, int py, MPI_Status* Status) {
	int neighbors[4];
	neighbors[2] = rank_y == 0 ? -1 : px * (rank_y - 1) + rank_x;
	neighbors[1] = rank_x == 0 ? -1 : px * rank_y + rank_x - 1;
	neighbors[0] = rank_y == py - 1 ? -1 : px * (rank_y + 1) + rank_x;
	neighbors[3] = rank_x == px - 1 ? -1 : px * rank_y + rank_x + 1;

	// std::string message = "Rank: " + std::to_string(px * rank_y + rank_x) + " ";
	// for (int i = 0; i < 4; i++) {
	// 	message += std::to_string(neighbors[i]) + " ";
	// }
	// std::cout << message << std::endl;

	std::vector<double> buffer_x1_send((stop_x - start_x) * fict * 2);
	std::vector<double> buffer_x2_send((stop_x - start_x) * fict * 2);
	std::vector<double> buffer_y1_send((stop_y - start_y) * fict * 2);
	std::vector<double> buffer_y2_send((stop_y - start_y) * fict * 2);

	std::vector<double> buffer_x1_recv((stop_x - start_x) * fict * 2);
	std::vector<double> buffer_x2_recv((stop_x - start_x) * fict * 2);
	std::vector<double> buffer_y1_recv((stop_y - start_y) * fict * 2);
	std::vector<double> buffer_y2_recv((stop_y - start_y) * fict * 2);

	if (neighbors[0] != -1) {
		for (int i = stop_y - start_y; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				buffer_x1_send[(i - stop_y + start_y) * (stop_x - start_x) + j - fict] = u[i][j];
				buffer_x1_send[(stop_x - start_x) * fict + ((i - stop_y + start_y) * (stop_x - start_x) + j - fict)] = v[i][j];
			}
		}
	}
	else {
		for (int i = stop_y - start_y; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
				v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
			}
		}
	}

	// std::cout << "first send" << std::endl;
	
	if (neighbors[1] != -1) {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < fict + fict; j++) {
				buffer_y1_send[(i - fict) * fict + j - fict] = u[i][j];
				buffer_y1_send[(stop_y - start_y) * fict + ((i - fict) * fict + j - fict)] = v[i][j];
			}	
		}
	}
	else {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < fict + fict; j++) {
				u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
				v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
			}	
		}
	}


	// std::cout << "second send" << std::endl;
	
	if (neighbors[2] != -1) {
		for (int i = fict; i < fict + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				buffer_x2_send[(i - fict) * (stop_x - start_x) + j - fict] = u[i][j];
				buffer_x2_send[(stop_x - start_x) * fict + ((i - fict) * (stop_x - start_x) + j - fict)] = v[i][j];
			}
		}
	}
	else {
		for (int i = fict; i < fict + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
				v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
			}
		}
	}

	// std::cout << "third send" << std::endl;

	if (neighbors[3] != -1) {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = stop_x - start_x; j < stop_x - start_x + fict; j++) {
				buffer_y2_send[(i - fict) * fict + j - stop_x + start_x] = u[i][j];
				buffer_y2_send[(stop_y - start_y) * fict + ((i - fict) * fict + j - stop_x + start_x)] = v[i][j];
			}
		}
	}
	else {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = stop_x - start_x; j < stop_x - start_x + fict; j++) {
				u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
				v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., time, Re);
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank_x != px - 1) {
		MPI_Send(buffer_y2_send.data(), (stop_y - start_y) * fict * 2, MPI_DOUBLE, neighbors[3], 3, MPI_COMM_WORLD); // right
	}
	if (rank_x != 0) {
		MPI_Recv(buffer_y1_recv.data(), (stop_y - start_y) * fict * 2, MPI_DOUBLE, neighbors[1], 3, MPI_COMM_WORLD, Status); // left
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank_x != 0) {
		MPI_Send(buffer_y1_send.data(), (stop_y - start_y) * fict * 2, MPI_DOUBLE, neighbors[1], 1, MPI_COMM_WORLD); // left
	}
	if (rank_x != px - 1) {
		MPI_Recv(buffer_y2_recv.data(), (stop_y - start_y) * fict * 2, MPI_DOUBLE, neighbors[3], 1, MPI_COMM_WORLD, Status); // right
	}
	MPI_Barrier(MPI_COMM_WORLD);


	if (rank_y != py - 1) {
		MPI_Send(buffer_x1_send.data(), (stop_x - start_x) * fict * 2, MPI_DOUBLE, neighbors[0], 0, MPI_COMM_WORLD); // up
	}
	if (rank_y != 0) {
		MPI_Recv(buffer_x2_recv.data(), (stop_x - start_x) * fict * 2, MPI_DOUBLE, neighbors[2], 0, MPI_COMM_WORLD, Status); // down
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank_y != 0) {
		MPI_Send(buffer_x2_send.data(), (stop_x - start_x) * fict * 2, MPI_DOUBLE, neighbors[2], 2, MPI_COMM_WORLD); // down
	}
	if (rank_y != py - 1) {
		MPI_Recv(buffer_x1_recv.data(), (stop_x - start_x) * fict * 2, MPI_DOUBLE, neighbors[0], 2, MPI_COMM_WORLD, Status); // up
	}
	MPI_Barrier(MPI_COMM_WORLD);


	
	// std::cout << "fourth send" << std::endl;
	
	if (neighbors[0] != -1) {
		for (int i = stop_y - start_y + fict; i < stop_y - start_y + fict + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				u[i][j] = buffer_x1_recv[(i - stop_y + start_y - fict) * (stop_x - start_x) + j - fict];
				v[i][j] = buffer_x1_recv[(stop_x - start_x) * fict + ((i - stop_y + start_y - fict) * (stop_x - start_x) + j - fict)];
			}
		}
	}

	
	// std::cout << "first recv" << std::endl;
	
	if (neighbors[1] != -1) {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = 0; j < fict; j++) {
				u[i][j] = buffer_y1_recv[(i - fict) * fict + j];
				v[i][j] = buffer_y1_recv[(stop_y - start_y) * fict + ((i - fict) * fict + j)];
			}	
		}
	}

	// std::cout << "second recv" << std::endl;
	
	if (neighbors[2] != -1) {
		for (int i = 0; i < fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				u[i][j] = buffer_x2_recv[i * (stop_x - start_x) + j - fict];
				v[i][j] = buffer_x2_recv[(stop_x - start_x) * fict + (i * (stop_x - start_x) + j - fict)];
			}
		}
	}

	// std::cout << "third recv" << std::endl;

	if (neighbors[3] != -1) {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = stop_x - start_x + fict; j < stop_x - start_x + fict + fict; j++) {
				u[i][j] = buffer_y2_recv[(i - fict) * fict + j - stop_x + start_x - fict];
				v[i][j] = buffer_y2_recv[(stop_y - start_y) * fict + ((i - fict) * fict + j - stop_x + start_x - fict)];
			}
		}
	}

	// std::cout << "fourth recv" << std::endl;
}


void write_data_parallel(std::vector<std::vector<double>> u, std::vector<std::vector<double>> v, double time, double Re, double dx, double dy,
					int start_x, int stop_x, int start_y, int stop_y, 
					int fict, int rank_x, int rank_y, std::string filename) {
    if (rank_x + rank_y == 0) {
        if (!fs::exists(filename)) {
            fs::create_directory(filename);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    std::string outpath = filename + "/" + "Time=" + std::to_string(time) + ".csv";
    std :: ostringstream buffer;
    
    if (rank_x + rank_y == 0) {
        buffer << "X,Y,U_exp,V_exp,U_theory,V_theory\n";
    }
	double x, y;
    for (int i = fict; i < stop_y - start_y + fict; i++) {
        for (int j = fict; j < stop_x - start_x + fict; j++) {
			x = (start_x + (j - fict)) * dx + dx / 2.;
			y = (start_y + (i - fict)) * dy + dy / 2.;
        	buffer << x << "," << y << "," << u[i][j] << "," << v[i][j] << "," << u_theory(x, y, time, Re) << "," << v_theory(x, y, time, Re) << "\n";
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
	double t_end = 3. / ddt / ddt;
	int Nx = N;
	int Ny = N;
	double dt = 0.001 / ddt / ddt;
	double Re = 100;
	int fict = 1;
	double dx = 1. / Nx;
	double dy = 1. / Ny;
	std::string message;
	std::string filename = "Out_csvs_N_" + std::to_string(N);
	int iterwrite = 30;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status Status;

	filename += "_p_" + std::to_string(size);

	int px, py;
	for (int i = sqrt(size); i > 0; i--) {
		if (size % i == 0) {
			px = i;
			py = size / i;
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

	// message = "Rank " + std::to_string(rank) + " start_x: " + std::to_string(start_x) + " stop_x: " + std::to_string(stop_x) + " start_y: " + std::to_string(start_y) + + " stop_y: " + std::to_string(stop_y);
	// std::cout << message << std::endl;

	double manytime = MPI_Wtime();
	std::vector<std::vector<double>> u(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	std::vector<std::vector<double>> v(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	std::vector<std::vector<double>> new_u(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	std::vector<std::vector<double>> new_v(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	
	InitialState(u, v, start_x, stop_x, start_y, stop_y, fict, dx, dy, Re);
	double tmp_Re_part = 0;

	// message = "Rank " + std::to_string(rank) + " done initializing";
	// std::cout << message << std::endl;

	int iter = 1;
	write_data_parallel(u, v, t_start, Re, dx, dy, start_x, stop_x, start_y, stop_y, fict, rank_x, rank_y, filename);
	while (t_start < t_end) {
		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				tmp_Re_part = ((u[i][j + 1] - 2 * u[i][j] + u[i][j - 1]) / dx / dx + (u[i + 1][j] - 2. * u[i][j] + u[i - 1][j]) / dy / dy) / Re;
				new_u[i][j] = u[i][j] + (tmp_Re_part - u[i][j] * (u[i][j + 1] - u[i][j - 1]) / 2. / dx - v[i][j] * (u[i + 1][j] - u[i - 1][j]) / 2. / dy) * dt;
				tmp_Re_part = ((v[i][j + 1] - 2 * v[i][j] + v[i][j - 1]) / dx / dx + (v[i + 1][j] - 2. * v[i][j] + v[i - 1][j]) / dy / dy) / Re;
				new_v[i][j] = v[i][j] + (tmp_Re_part - u[i][j] * (v[i][j + 1] - v[i][j - 1]) / 2. / dx - v[i][j] * (v[i + 1][j] - v[i - 1][j]) / 2. / dy) * dt;
			}
		}

		for (int i = fict; i < stop_y - start_y + fict; i++) {
			for (int j = fict; j < stop_x - start_x + fict; j++) {
				u[i][j] = new_u[i][j];
				v[i][j] = new_v[i][j];
			}
		}

		t_start += dt;
		// message = "Rank_x " + std::to_string(rank_x) + " rank_y " + std::to_string(rank_y) + " reached boundary";
		// std::cout << message << std::endl;
		Boundary(u, v, t_start, Re, dx, dy, start_x, stop_x, start_y, stop_y, fict, rank_x, rank_y, px, py, &Status);
		MPI_Barrier(MPI_COMM_WORLD);

		//message = "";
		//message += "I am proc" + std::to_string(rank) + " time: " + std::to_string(t_start);
		//message += " x_coord: " + std::to_string(start_x * dx + dx / 2.) + " y_coord: " + std::to_string(start_y * dx + dx / 2.);
		//message += " u_theory: " + std::to_string(u_theory(start_x * dx + dx / 2., start_y * dx + dx / 2., t_start, Re)) + " v_theory: " + std::to_string(v_theory(start_x * dx + dx / 2., start_y * dx + dx / 2., t_start, Re));
		//message += " u_exp: " + std::to_string(u[fict][fict]) + " v_exp: " + std::to_string(v[fict][fict]);
		//std::cout << message << std::endl;

		if (iter % iterwrite == 0 || t_start == t_end) {
			write_data_parallel(u, v, t_start, Re, dx, dy, start_x, stop_x, start_y, stop_y, fict, rank_x, rank_y, filename);
		}
		iter++;
	}

	// message = "";
	// message += "I am proc" + std::to_string(rank);
	// message += " x_coord: " + std::to_string(start_x * dx + dx / 2.) + " y_coord: " + std::to_string(start_y * dx + dx / 2.);
	// message += " u_theory: " + std::to_string(u_theory(start_x * dx + dx / 2., start_y * dx + dx / 2., t_start, Re)) + " v_theory: " + std::to_string(v_theory(start_x * dx + dx / 2., start_y * dx + dx / 2., t_start, Re));
	// message += " u_exp: " + std::to_string(u[fict][fict]) + " v_exp: " + std::to_string(v[fict][fict]);
	// std::cout << message << std::endl;
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		std::cout << "\nTime of working: " << MPI_Wtime() - manytime << "\n";
	}
	MPI_Finalize();
	return 0;
}
