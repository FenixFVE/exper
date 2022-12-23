#pragma once

#include "grid.hpp"
#include "sparse_matrix.hpp"

void grid_test() {
	Grid grid("grid.txt");

	cout << "X: ";
	for (auto x : grid.grid_x) cout << x << " ";
	cout << '\n';

	cout << "Y: ";
	for (auto y : grid.grid_y) cout << y << " ";
	cout << '\n';

	cout << "Number to coordinates and back:\n";
	for (int i = 0; i < grid.number_to_coordinates.size(); ++i) {
		auto [x, y] = grid.number_to_coordinates[i];
		cout << i << ": " << x << ", " << y << " | " << grid.coordinates_to_number[{ x, y }] << "\n";
	}

	cout << "Elements: \n";
	for (int i = 0; i < grid.elements.size(); ++i) {
		cout << i << ": ";
		for (auto x : grid.elements[i]) {
			cout << x << " ";
		}
		cout << "\n";
	}

	cout << "Adj: \n";
	for (int i = 0; i < grid.adjacency_list.size(); ++i) {
		cout << i << ": ";
		for (int x : grid.adjacency_list[i]) cout << x << " ";
		cout << '\n';
	}
}

void MSG_test() {
	SparseMatrix matrix;
	matrix.n = 5;
	matrix.di = { 1,1,1,1,1 };
	matrix.jg = { 0,0,1,2,2,3 };
	matrix.ig = { 0,0,1,3,4,6 };
	matrix.gg = { 1,2,3,4,5,6 };
	vector<double> b = { 4,5,15,11,12 };
	vector<double> x0 = { 0,0,0,0,0 };
	auto x = MSG_solve(matrix, x0, b, 1e-16, 1000);
	for (auto i : x) cout << i << '\n';
}