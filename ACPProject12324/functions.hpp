#include <vector>
#include <cmath>
#include <random>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <execution>
#include <mutex>
#include <thread>
#include <boost/unordered_map.hpp>

double pi = atan(1) * 4;

double positive_mod(double a, double b) {
	return std::fmod((b + std::fmod(a, b)), b);
}
//distance function
double* min_image_distance(double x1, double y1, double x2, double y2, const double L) {

	if (x1 > x2 + L / 2) {
		x2 += L;
	}
	else if (x1 < x2 - L / 2) {
		x2 -= L;
	}
	if (y1 > y2 + L / 2) {
		y2 += L;
	}
	else if (y1 < y2 - L / 2) {
		y2 -= L;
	}
	double d[2]{};
	d[0] = x1 - x2;
	d[1] = y1 - y2;

	return d;
}

//function to create initial configurations
std::vector<std::vector<double>> init_uniform_random_grid(const int N, const double rho, bool orientation) {

	std::vector<std::vector<double>> grid(2);
	std::vector<double> null(N);
	std::fill(null.begin(), null.end(), 0);
	std::fill(grid.begin(), grid.end(), null);

	double side_length = sqrt(N / rho);

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0, side_length);

	if (orientation) {
		grid.push_back(null);
		std::uniform_real_distribution<double> dist_phi(0, 2 * pi);

		for (int i = 0; i < N; i++) {
			grid[0][i] = dist(gen);
			grid[1][i] = dist(gen);
			grid[2][i] = dist_phi(gen);
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			grid[0][i] = dist(gen);
			grid[1][i] = dist(gen);
		}
	}
	return grid;
}
std::vector<std::vector<double>> init_uniform_grid(const int N, const double rho, bool orientation) {

	std::vector<std::vector<double>> grid(2);
	std::vector<double> null(N);
	std::fill(null.begin(), null.end(), 0);
	std::fill(grid.begin(), grid.end(), null);

	double side_length = sqrt(N / rho);
	int n = ceil(sqrt(N));
	double a = side_length / n;

	std::random_device rd;
	std::mt19937 gen(rd());

	if (orientation) {
		grid.push_back(null);
		std::uniform_real_distribution<double> dist_phi(0, 2 * pi);

		for (int i = 0; i < N; i++) {
			grid[0][i] = i % n * a;
			grid[1][i] = floor(i / n) * a;
			grid[2][i] = dist_phi(gen);
		}
	}
	else {
		for (int i = 0; i < N; i++) {
			grid[0][i] = i % n * a;
			grid[1][i] = floor(i / n) * a;
		}
	}
	return grid;
}



//algorithms to find pairs


//naive approach
std::vector<std::pair<int, int>> neighbors_naive(const std::vector<double>& x_values, const std::vector<double>& y_values, double d_max) {
	std::vector<std::pair<int, int>> neighbors;

	for (int i = 0; i < x_values.size() - 1; i++) {
		for (int j = i + 1; j < x_values.size(); j++) {
			double d2 = (x_values[i] - x_values[j]) * (x_values[i] - x_values[j]) + (y_values[i] - y_values[j]) * (y_values[i] - y_values[j]);
			if (d2 < d_max*d_max) {
				neighbors.push_back({ i, j });
			}
		}
	}
	return neighbors;
}

//cell list 
std::vector<std::pair<int, int>> cell_list_algo(std::vector<double>& x_values, std::vector<double>& y_values, const double d2_max, const double L, const int N, const int n, const double cell_length) {

	std::vector<std::pair<int, int>> pairs;

	// Create a hash map for cell-wise particle storage

	std::unordered_map<int, std::vector<int>> cell_map;
	cell_map.reserve(n * n);

	for (int i = 0; i < N; i++) {
		//cell_map[i].reserve(10); // reserving some memory for performance
		int cell_x = floor(x_values[i] / cell_length);
		int cell_y = floor(y_values[i] / cell_length);

		int cell_label = cell_x + cell_y * n; //finding cell label of the particle

		if (i == 0) {
			cell_map[cell_label].push_back(i); // first particle so dont need to check for pairs
		}
		else { //if not first particle checking neighboring cells for pairs
			for (int dx = -1; dx < 2; dx++) {
				for (int dy = -1; dy < 2; dy++) {
					int x_neighbor = cell_x + dx;
					int y_neighbor = cell_y + dy;

					if (0 <= x_neighbor && x_neighbor < n && 0 <= y_neighbor && y_neighbor < n) {
						int neighbor_label = x_neighbor + y_neighbor * n;
						for (int p2 : cell_map[neighbor_label]) { //checking each particle in neighbor cell

							double d2 = (x_values[i] - x_values[p2]) * (x_values[i] - x_values[p2]) + (y_values[i] - y_values[p2]) * (y_values[i] - y_values[p2]); //distance

							if (d2 < d2_max) { //if distance small enough add to the force
								pairs.push_back({ i,p2 });
							}
						}
					}
				}
			}
			cell_map[cell_label].push_back(i); //grouping particle in its cell
		}
	}
	return pairs;
}



//integrator for langevin equation without interactions
void langevin_integrator(const std::vector<double>& x_values, const std::vector<double>& y_values, const std::vector<double>& phi_values, const double L, const double t_end, const double pe, const char* path_msd, const char* path_trajectories) {

	std::ofstream ex2_msd(path_msd);
	std::ofstream ex2_trajectories(path_trajectories);

	int N = x_values.size();

	std::vector<std::vector<double>> initial_positions;
	std::vector<std::vector<double>> positions;
	std::vector<std::vector<double>> total_displacement(3);

	std::vector<double> null(N);
	std::fill(null.begin(), null.end(), 0);
	std::fill(total_displacement.begin(), total_displacement.end(), null);


	initial_positions = { x_values,y_values,phi_values };
	positions = initial_positions;

	int saved_particles[5] = { 100,300,500,700,900 };

	for (int i = 0; i < std::size(saved_particles); i++) {

		ex2_trajectories << "Particle(" << saved_particles[i] << ")   x   y   phi   ";

	}
	ex2_trajectories << "\n";

	for (int i = 0; i < std::size(saved_particles); i++) {
		ex2_trajectories << positions[0][saved_particles[i]] << " " << positions[1][saved_particles[i]] << " " << positions[2][saved_particles[i]] << " ";
	}
	ex2_trajectories << "\n";

	const double delta_t = 0.0001;
	const double q = 1.543;
	double xi_x;
	double xi_y;
	double chi;
	double dx;
	double dy;
	double dphi;
	double msd_value = 0;
	double msad_value = 0;
	std::minstd_rand gen(76);
	std::normal_distribution<double> dist(0, 1);

	for (int t = 0; t * delta_t < t_end; t++) {
		msd_value = 0;
		msad_value = 0;

		for (int i = 0; i < x_values.size(); i++) {

			xi_x = dist(gen);
			xi_y = dist(gen);
			chi = dist(gen);

			dx = pe * cos(positions[2][i]) * delta_t + sqrt(2 * delta_t) * xi_x;
			dy = pe * sin(positions[2][i]) * delta_t + sqrt(2 * delta_t) * xi_y;
			dphi = q * sqrt(2 * delta_t) * chi;

			positions[0][i] = positive_mod(positions[0][i] + dx, L);
			positions[1][i] = positive_mod(positions[1][i] + dy, L);
			positions[2][i] = positive_mod(positions[2][i] + dphi, 2 * pi);

			total_displacement[0][i] += dx;
			total_displacement[1][i] += dy;
			total_displacement[2][i] += dphi;

			msd_value += total_displacement[0][i] * total_displacement[0][i] + total_displacement[1][i] * total_displacement[1][i];
			msad_value += total_displacement[2][i] * total_displacement[2][i];

		}
		for (int j = 0; j < std::size(saved_particles); j++) {
			ex2_trajectories << positions[0][saved_particles[j]] << " " << positions[1][saved_particles[j]] << " " << positions[2][saved_particles[j]] << " ";
		}
		ex2_trajectories << "\n";
		ex2_msd << t * delta_t << " " << msd_value / N << " " << msad_value / N << "\n";
	}
}


//function that gets and stores the index of neighboring cells for cell list algo with periodic boundary conditions

std::unordered_map<int, std::array<int, 9>> cell_neighbors_periodic(const int n) {
	std::unordered_map<int, std::array<int, 9>> map_neighbors;
	map_neighbors.reserve(n * n);

	for (int i = 0; i < n * n; i++) {
		int x_label = i % n;
		int y_label = floor(i / n);
		int index = 0;

		for (int dx = -1; dx < 2; dx++) {
			for (int dy = -1; dy < 2; dy++) {
				int x_neighbor = (n + ((x_label + dx) % n)) % n;
				int y_neighbor = (n + ((y_label + dy) % n)) % n;
				int neighbor_label = x_neighbor + y_neighbor * n;
				map_neighbors[i][index] = neighbor_label;
				index++;
			}
		}
	}
	return map_neighbors;
}



// function for cell list algorithm that directly calculates and returns the resulting force on each particle 
std::vector<std::pair<double, double>> cell_force(std::vector<double>& x_values, std::vector<double>& y_values, std::unordered_map<int, std::array<int, 9>>& map_neighbors, const double d2_max, const double L,const int N,const int n ,const double cell_length) {

	std::pair<double, double> null = { 0,0 };
	std::vector<std::pair<double, double>> force(N);
	std::fill(force.begin(), force.end(), null);

	double f{};

	// Create a hash map for cell-wise particle storage
	std::unordered_map<int, std::vector<int>> cell_map;
	cell_map.reserve(n * n);

	for (int i = 0; i < N; i++) {
		uint16_t cell_x = floor(x_values[i] / cell_length);
		uint16_t cell_y = floor(y_values[i] / cell_length);

		uint16_t cell_label = cell_x + cell_y * n; //finding cell label of the particle

		if (i == 0) {
			cell_map[cell_label].push_back(i); // first particle so dont need to check for pairs
		}
		else { //if not first particle checking neighboring cells for pairs
			for (uint16_t neighbor_label : map_neighbors[cell_label]){
				for (uint16_t p2 : cell_map[neighbor_label]) { //checking each particle in neighbor cell
						
					double* distance = min_image_distance(x_values[i], y_values[i], x_values[p2], y_values[p2], L); //distance
					
					double dx = distance[0];
					double dy = distance[1];
					double d2 = dx * dx + dy * dy;
					if (d2 < d2_max) { //if distance small enough add to the force
						f = 1 / (d2 * d2 * d2 * d2 * d2 * d2 * d2) - 0.5 / (d2 * d2 * d2 * d2);
						force[i].first += f * dx;
						force[i].second += f * dy;
						force[p2].first -= f * dx;
						force[p2].second -= f * dy;
					}	
				}
			}
			cell_map[cell_label].push_back(i); //grouping particle in its cell
		}	
	}
	return force;
}

//langevin integrator with interactions  

void integrator_interactions(const std::vector<double>& x_values, const std::vector<double>& y_values, const std::vector<double>& phi_values, const double L, const double t_end, const double pe, const double d_max, double d2_max, const char* path_msd,const char* snapshots) {

	std::ofstream off_msd(path_msd);
	std::ofstream off_snap(snapshots);

	const double delta_t = 0.00001;
	int nt = t_end / delta_t + 1;
	const double q = 1.543;
	const int N = x_values.size();
	double d_x;
	double d_y;
	double d_phi;
	double xi_x;
	double xi_y;
	double chi;
	double msd_value = 0;
	double msad_value = 0;
	std::vector<double> x_list = x_values;
	std::vector<double> y_list = y_values;
	std::vector<double> phi_list = phi_values;
	std::vector<double> x_displacement(N); 
	std::vector<double> y_displacement(N);
	std::vector<double> phi_displacement(N);
	std::vector<double> null(N);
	std::fill(null.begin(), null.end(), 0);
	x_displacement = null;
	y_displacement = null;
	phi_displacement = null;

	std::minstd_rand gen(76);
	std::normal_distribution<double> dist(0, 1);

	const int n = floor(L / d_max);
	double cell_length = L / n;
	std::unordered_map<int, std::array<int, 9>> map_neighbors = cell_neighbors_periodic(n);


	//snapshot at t=0
	off_snap << 0 << " ";
	for (int i = 0; i < N; i++) {
		off_snap << x_list[i] << " " << y_list[i] << " " << phi_list[i] << " ";
	}
	off_snap << "\n";

	for (int t = 1; t <= nt; t++) {
		msd_value = 0;
		msad_value = 0;

		
		std::vector<std::pair<double, double>> f = cell_force(x_list, y_list,map_neighbors, d2_max, L,N,n,cell_length);

		for (int i = 0; i < N; i++) {
			xi_x = dist(gen);
			xi_y = dist(gen);
			chi = dist(gen);

			d_phi = q * sqrt(2 * delta_t) * chi;
			d_x = pe * cos(phi_list[i]) * delta_t + sqrt(2 * delta_t) * xi_x + delta_t * delta_t * f[i].first;
			d_y = pe * sin(phi_list[i]) * delta_t + sqrt(2 * delta_t) * xi_y + delta_t * delta_t * f[i].second;

			phi_list[i] = positive_mod(phi_list[i] + d_phi, 2 * pi);
			x_list[i] = positive_mod(x_list[i] + d_x, L);
			y_list[i] = positive_mod(y_list[i] + d_y, L);

			x_displacement[i] += d_x;
			y_displacement[i] += d_x;
			phi_displacement[i] += d_phi;

			msd_value += x_displacement[i] * x_displacement[i] + y_displacement[i] * y_displacement[i];
			msad_value += phi_displacement[i] * phi_displacement[i];
		}
		if ( (t) % 10000 == 0) {
			off_snap << (t)*delta_t << " ";
			for (int i = 0; i < N; i++) {
				off_snap << x_list[i] << " " << y_list[i] << " " << phi_list[i] << " ";
			}
			off_snap << "\n";
		}

		// here output the msd
		if (t % 1000 == 0) {
			off_msd << t * delta_t << " " << msd_value / N << " " << msad_value / N << "\n";
		}	
	}
}