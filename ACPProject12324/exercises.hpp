#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <chrono>
#include <algorithm>
#include <execution>
//#include "functions.hpp"


//Exercise 1 a)

void exercise1a() {
	//define example values 1a
	std::ofstream ex1a_output("ex1a.dat");
	std::vector<double> x_values = { 0,1,1,2 };
	std::vector<double> y_values = { 0,1,0,1 };
	double d_max = 1.5;

	std::vector<std::pair<int, int>> pairs = neighbors_naive(x_values, y_values, d_max);

	for (auto i : pairs) {
		ex1a_output << i.first << " " << i.second << std::endl;
	}
}

//Exercise 1 b)
void exercise1b() {

	std::ofstream grid_output("ex1b_grid.dat");
	std::ofstream pairs_naive15_output("ex1b_pairs_naive15.dat");
	std::ofstream pairs_naive13_output("ex1b_pairs_naive13.dat");


	//creating the 25x25 grid 
	std::vector<double> x_values_b;
	std::vector<double> y_values_b;

	for (int i = 0; i < 25; i++) {
		for (int j = 0; j < 25; j++) {
			x_values_b.push_back(i);
			y_values_b.push_back(j);
		}
	}
	//finding the pairs for the two different d_max
	double d_max_b15 = 1.5;
	double d_max_b13 = 1.3;

	std::vector<std::pair<int, int>> pairs_b15 = neighbors_naive(x_values_b, y_values_b, d_max_b15);
	std::vector<std::pair<int, int>> pairs_b13 = neighbors_naive(x_values_b, y_values_b, d_max_b13);


	//output of the grid and pairs
	for (int i = 0; i < x_values_b.size(); i++) {
		grid_output << x_values_b[i] << " " << y_values_b[i] << std::endl;
	}
	for (int i = 0; i < pairs_b15.size(); i++) {
		pairs_naive15_output << pairs_b15[i].first << " " << pairs_b15[i].second << std::endl;
	}
	for (int i = 0; i < pairs_b13.size(); i++) {
		pairs_naive13_output << pairs_b13[i].first << " " << pairs_b13[i].second << std::endl;
	}

	std::cout << "for 1.5 n=" << pairs_b15.size() << " and for 1.3 n=" << ' ' << pairs_b13.size() << std::endl;
	// for 1.5 dmax counting the connections goes like 24*(3+23*4+2)+24
	// for 1.3 dmax counting the connections goes like 24*(2*24+1)+24
}

//exercise 1 c) (and 1e using the cell list algorithm on the same grid as the naive algorithm)

void exercise1c() {
	std::ofstream ex1c_grid("ex1c_grid.dat");
	std::ofstream ex1c_pairs_naive("ex1c_pairs_naive.dat");
	std::ofstream ex1c_pairs_cell("ex1c_pairs_cell.dat");

	int N = 500;
	double L = 1;
	double rho = N/L; // side length is 1 and N=500 so number_density rho=500
	bool orientation = false;
	std::vector<std::vector<double>> x_y_values = init_uniform_random_grid(N, rho, orientation);

	double d_max = 0.08;
	const int n = floor(L / d_max);
	const double cell_length = L / n;

	std::vector<std::pair<int, int>> pairs_naive = neighbors_naive(x_y_values[0], x_y_values[1], d_max);

	std::vector<std::pair<int, int>> pairs_cell = cell_list_algo(x_y_values[0], x_y_values[1],d_max*d_max,L,N,n,cell_length);

	for (int i = 0; i < x_y_values[0].size(); i++) {
		ex1c_grid << x_y_values[0][i] << " " << x_y_values[1][i] << std::endl;
	}
	for (int i = 0; i < pairs_naive.size(); i++) {
		ex1c_pairs_naive << pairs_naive[i].first << " " << pairs_naive[i].second << std::endl;
	}
	for (int i = 0; i < pairs_cell.size(); i++) {
		ex1c_pairs_cell << pairs_cell[i].first << " " << pairs_cell[i].second << std::endl;
	}
}

//exercise 1 e) (redo part 1b with the cell list algo)

void exercise1e() {
	std::ofstream pairs_cell15_output("ex1e_pairs_cell15.dat");
	std::ofstream pairs_cell13_output("ex1e_pairs_cell13.dat");


	//creating the 25x25 grid 
	std::vector<double> x_values;
	std::vector<double> y_values;

	for (int i = 0; i < 25; i++) {
		for (int j = 0; j < 25; j++) {
			x_values.push_back(i);
			y_values.push_back(j);
		}
	}
	//finding the pairs for the two different d_max
	double d_max_15 = 1.5;
	double d_max_13 = 1.3;

	std::vector<std::pair<int, int>> pairs_15 = neighbors_naive(x_values, y_values, d_max_15);
	std::vector<std::pair<int, int>> pairs_13 = neighbors_naive(x_values, y_values, d_max_13);


	//output of the pairs
	for (int i = 0; i < pairs_15.size(); i++) {
		pairs_cell15_output << pairs_15[i].first << " " << pairs_15[i].second << std::endl;
	}
	for (int i = 0; i < pairs_13.size(); i++) {
		pairs_cell13_output << pairs_13[i].first << " " << pairs_13[i].second << std::endl;
	}

	std::cout << "for 1.5 n=" << pairs_15.size() << " and for 1.3 n=" << ' ' << pairs_13.size() << std::endl;
	// for 1.5 dmax counting the connections goes like 24*(3+23*4+2)+24
	// for 1.3 dmax counting the connections goes like 24*(2*24+1)+24
}

//exercise 1 f) 

void exercise1f() {

	std::ofstream ex1f_runtime_output("ex1f_runtime.dat");

	std::vector<int> system_sizes = { 10 , 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000,25000,50000,75000,100000 };
	
	double side_length;
	double d_max = 1.5;
	double d2_max = d_max * d_max;
	bool orientation = false;
	std::vector<std::vector<double>> x_y_values;
	


	for (int i : system_sizes) {
		side_length = sqrt(i);
		const int n = floor(side_length / d_max);
		const double cell_length = side_length / n;
		x_y_values = init_uniform_grid(i, 1, orientation);


		auto start = std::chrono::high_resolution_clock::now();

		std::vector<std::pair<int, int>> pairs_naive = neighbors_naive(x_y_values[0], x_y_values[1], d_max);

		auto stop = std::chrono::high_resolution_clock::now();

		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);



		auto start1 = std::chrono::high_resolution_clock::now();

		std::vector<std::pair<int, int>> pairs_cell = cell_list_algo(x_y_values[0], x_y_values[1],d2_max,side_length,i,n,cell_length);
		
		auto stop1 = std::chrono::high_resolution_clock::now();

		auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);

		ex1f_runtime_output << i << " " << duration.count() << " " << duration1.count() << "\n";	
	}
}


//exercise 2 

void exercise2() {

	const char* offstream_traj[2] = { "ex2_trajectories_0.dat" ,"ex2_trajectories_20.dat"  };
	const char* offstream_msd[2] = { "ex2_msd_0.dat" ,"ex2_msd_20.dat" };

	int N = 1000;
	double L = 85;
	double rho = N / (L * L);
	bool orientation = true;
	double t_end = 100;
	double Pe[2] = { 0,20 };
	int iterator[2] = { 0,1 };

	std::vector<std::vector<double>> grid = init_uniform_grid(N, rho, orientation);

	std::for_each(std::execution::par, std::begin(iterator), std::end(iterator), [grid, L, t_end, offstream_traj,offstream_msd, Pe](int i) {

		langevin_integrator(grid[0], grid[1], grid[2], L, t_end, Pe[i], offstream_msd[i], offstream_traj[i]);
	});
}

//exercise 3

void exercise3() {
	int N = 1000;
	double L = 85;
	double rho = N / (L * L);
	bool orientation = true;
	double t_end = 100;
	double d_max = 2.5;
	double r2_cut = 2.5 * 2.5;

	double Pe[2] = { 0,20 };
	int iterator[2] = { 0,1 };


	const char* paths[2] = { "ex3_msd_0.dat" , "ex3_msd_20.dat" };
	const char* snapshots[2] = { "ex3_snap_0.dat", "ex3_snap_20.dat" };

	std::vector<std::vector<double>> grid = init_uniform_grid(N, rho, orientation);
	
	std::for_each(std::execution::par,std::begin(iterator), std::end(iterator), [grid,L,t_end,d_max,r2_cut,paths,Pe,snapshots](int i) {
		
		integrator_interactions(grid[0], grid[1], grid[2], L, t_end, Pe[i], d_max, r2_cut, paths[i],snapshots[i]);
	});
}
void exercise4() {
	int N = 4000;
	double L = 85;
	double rho = N / (L * L);
	bool orientation = true;
	double t_end = 100;
	double d_max = 2.5; 
	double d2_max = d_max*d_max;

	double Pe[6] = { 0,20,50,80 ,130,180 };

	int iterator[6] = { 0,1,2,3 ,4,5};


	const char* paths[6] = { "ex4_msd_0.dat" , "ex4_msd_20.dat", "ex4_msd_50.dat" , "ex4_msd_80.dat","ex4_msd_130.dat","ex4_msd_180.dat" };
	const char* snapshots[6] = { "ex4_snap_0.dat", "ex4_snap_20.dat", "ex4_snap_50.dat" ,"ex4_snap_80.dat","ex4_snap_130.dat","ex4_snap_180.dat" };

	std::vector<std::vector<double>> grid = init_uniform_grid(N, rho, orientation);

	std::for_each(std::execution::par, std::begin(iterator), std::end(iterator), [grid, L, t_end, d_max, d2_max, paths, Pe,snapshots](int i) {

		integrator_interactions(grid[0], grid[1], grid[2], L, t_end, Pe[i], d_max, d2_max, paths[i], snapshots[i]);
		});
}
