#pragma once

#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <map>
#include <iostream>

using namespace std;

class Grid {
public:
	// Readonly
	int x_size;
	int y_size;
	int grid_x_size;
	int grid_y_size;
	vector<double> grid_x;
	vector<double> grid_y;
	int nodes_number;
	vector<pair<int,int>> number_to_coordinates;
	map<pair<int, int>, int> coordinates_to_number;
	int elem_number;
	int connections;
	vector <vector<int>> elements;
	vector<set<int>> adjacency_list;


	Grid(string grid_file_name) {

		ifstream grid_file(grid_file_name);
		if (!grid_file.good()) throw "Grid file error";
		double temporal;
		grid_file >> x_size >> y_size;
		
		grid_x_size = x_size * 2 - 1;
		grid_x.reserve(grid_x_size);
		for (int i = 0; i < x_size; ++i) {
			grid_file >> temporal;
			if (i != 0) {
				grid_x.push_back((temporal - grid_x.back()) / 2.0 + grid_x.back());
			}
			grid_x.push_back(temporal);
		}
		
		grid_y_size = y_size * 2 - 1;
		grid_y.reserve(grid_y_size);
		for (int j = 0; j < y_size; ++j) {
			grid_file >> temporal;
			if (j != 0) {
				grid_y.push_back((temporal - grid_y.back()) / 2.0 + grid_y.back());
			}
			grid_y.push_back(temporal);
		}
		
		int coordinates_counter = 0;
		number_to_coordinates.reserve(grid_x_size * grid_y_size);
		for (int y = 0; y < grid_y_size; ++y) {
			for (int x = 0; x < grid_x_size; ++x) {
				number_to_coordinates.push_back({ x, y });
				coordinates_to_number.insert({ {x,y}, coordinates_counter++ });
			}
		}
		
		int elements_counter = 0;
		elem_number = (x_size - 1) * (y_size - 1);
		elements.reserve(elem_number);
		for (int y = 0; y < grid_y_size - 1; y += 2) {
			for (int x = 0; x < grid_x_size - 1; x += 2) {
				vector<int> numbers;
				numbers.reserve(9);
				for (int i = 0; i < 3; ++i) {
					for (int j = 0; j < 3; ++j) {
						numbers.push_back(coordinates_to_number[{ x + j, y + i }]);
					}
				}
				elements.push_back(numbers);
			}
		}

		adjacency_list.resize(grid_x_size * grid_y_size, set<int>{});
		for (int i = 0; i < elements.size(); ++i) {
			for (int x : elements[i]) {
				for (int y : elements[i]) {
					if (y < x) {
						adjacency_list[x].insert(y);
					}
				}
			}
		}

		connections = 0;
		for (auto adj : adjacency_list) {
			connections += (int)adj.size();
		}

		nodes_number = grid_x_size * grid_y_size;
	}
};