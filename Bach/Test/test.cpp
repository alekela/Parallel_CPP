#include <iostream>
// #include <mpi.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>


double u_theory(double x, double y, double t, double Re) {
	return 0.75 - 0.25 * (1. / (1 + exp((-4. * x + 4. * y - t) * Re / 32.)));
}


double v_theory(double x, double y, double t, double Re) {
	return 0.75 + 0.25 * (1. / (1 + exp((-4. * x + 4. * y - t) * Re / 32.)));
}


int main(int argc, char* argv[]) {
	double t_start = 0;
	double t_end = 3;
	int Nx = 100;
	int Ny = 100;
	double dt = 0.001;
	double Re = 100;
	int fict = 1;
	double dx = 1. / Nx;
	double dy = 1. / Ny;
	std::string message;
	std::string filename = "Out_csvs";

    int size = 3;
    
    for (int rank = 0; rank < size; rank++) {

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
        
        std::vector<std::vector<double>> u(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	    std::vector<std::vector<double>> v(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	    std::vector<std::vector<double>> new_u(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
	    std::vector<std::vector<double>> new_v(stop_y - start_y + 2 * fict, std::vector<double>(stop_x - start_x + 2 * fict));
        
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

        std::vector<double> buffer_x1((stop_x - start_x) * fict * 2);
        std::vector<double> buffer_x2((stop_x - start_x) * fict * 2);
        std::vector<double> buffer_y1((stop_y - start_y) * fict * 2);
        std::vector<double> buffer_y2((stop_y - start_y) * fict * 2);

        std::vector<double> buffer_x1_recv((stop_x - start_x) * fict * 2);
        std::vector<double> buffer_x2_recv((stop_x - start_x) * fict * 2);
        std::vector<double> buffer_y1_recv((stop_y - start_y) * fict * 2);
        std::vector<double> buffer_y2_recv((stop_y - start_y) * fict * 2);

        if (neighbors[0] != -1) {
            for (int i = stop_y - start_y; i < stop_y - start_y + fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    buffer_x1[(i - stop_y + start_y) * (stop_x - start_x) + j - fict] = u[i][j];
                    buffer_x1[2 * ((i - stop_y + start_y) * (stop_x - start_x) + j - fict)] = v[i][j];
                }
            }
        }
        else {
            for (int i = stop_y - start_y; i < stop_y - start_y + fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                    v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                }
            }
        }

        // std::cout << "first send" << std::endl;
        
        if (neighbors[1] != -1) {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = fict; j < fict + fict; j++) {
                    buffer_y1[(i - fict) * fict + j - fict] = u[i][j];
                    buffer_y1[2 * ((i - fict) * fict + j - fict)] = v[i][j];
                }	
            }
        }
        else {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = fict; j < fict + fict; j++) {
                    u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                    v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                }	
            }
        }


        // std::cout << "second send" << std::endl;
        
        if (neighbors[2] != -1) {
            for (int i = fict; i < fict + fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    buffer_x2[(i - fict) * (stop_x - start_x) + j - fict] = u[i][j];
                    buffer_x2[2 * ((i - fict) * (stop_x - start_x) + j - fict)] = v[i][j];
                }
            }
        }
        else {
            for (int i = fict; i < fict + fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                    v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                }
            }
        }

        // std::cout << "third send" << std::endl;

        if (neighbors[3] != -1) {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = stop_x - start_x; j < stop_x - start_x + fict; j++) {
                    buffer_y2[(i - fict) * fict + j - stop_x + start_x] = u[i][j];
                    buffer_y2[2 * ((i - fict) * fict + j - stop_x + start_x)] = v[i][j];
                }
            }
        }
        else {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = stop_x - start_x; j < stop_x - start_x + fict; j++) {
                    u[i][j] = u_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                    v[i][j] = v_theory((start_x + j - fict) * dx + dx / 2., (start_y + i - fict) * dy + dy / 2., 0, Re);
                }
            }
        }


        
        // std::cout << "fourth send" << std::endl;
        
        if (neighbors[0] != -1) {
            for (int i = stop_y - start_y + fict; i < stop_y - start_y + fict + fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    u[i][j] = buffer_x1_recv[(i - stop_y + start_y - fict) * (stop_x - start_x) + j - fict];
                    v[i][j] = buffer_x1_recv[2 * ((i - stop_y + start_y - fict) * (stop_x - start_x) + j - fict)];
                }
            }
        }

        
        // std::cout << "first recv" << std::endl;
        
        if (neighbors[1] != -1) {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = 0; j < fict; j++) {
                    u[i][j] = buffer_y1_recv[(i - fict) * fict + j];
                    v[i][j] = buffer_y1_recv[2 * ((i - fict) * fict + j)];
                }	
            }
        }


        // std::cout << "second recv" << std::endl;
        
        if (neighbors[2] != -1) {
            for (int i = 0; i < fict; i++) {
                for (int j = fict; j < stop_x - start_x + fict; j++) {
                    u[i][j] = buffer_x2_recv[i * (stop_x - start_x) + j - fict];
                    v[i][j] = buffer_x2_recv[2 * (i * (stop_x - start_x) + j - fict)];
                }
            }
        }



        // std::cout << "third recv" << std::endl;

        if (neighbors[3] != -1) {
            for (int i = fict; i < stop_y - start_y + fict; i++) {
                for (int j = stop_x - start_x + fict; j < stop_x - start_x + fict + fict; j++) {
                    u[i][j] = buffer_y2_recv[(i - fict) * fict + j - stop_x + start_x - fict];
                    v[i][j] = buffer_y2_recv[2 * ((i - fict) * fict + j - stop_x + start_x - fict)];
                }
            }
        }

        // std::cout << "fourth recv" << std::endl;
    }
	return 0;
}
