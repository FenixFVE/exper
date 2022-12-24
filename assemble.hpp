#pragma once

#include "grid.hpp"
#include "sparse_matrix.hpp"
#include "parameters.hpp"
#include <cmath>
#include <algorithm>

class Matrix {
public:
	int rows, columns;
	vector<vector<double>> T;

	Matrix() { 
		rows = 0; columns = 0;
		T.reserve(1);
	}
	Matrix(int n, int m) {
		rows = n; columns = m;
		T = vector<vector<double>>(n, vector<double>(m, 0.0));
	}
	Matrix(const vector<vector<double>>& prot) {
		rows = prot.size(); columns = prot[0].size();
		T.resize(rows, vector<double>(columns, 0.0));
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				T[i][j] = prot[i][j];
			}
		}
	}
};

Matrix operator+(const Matrix& m1, const Matrix& m2) {
	if (m1.rows != m2.rows || m1.columns != m2.columns) {
		throw "matrix+ error";
	}
	Matrix n(m1.rows, m1.columns);
	for (int i = 0; i < m1.rows; ++i) {
		for (int j = 0; j < m1.columns; ++j) {
			n.T[i][j] = m1.T[i][j] + m2.T[i][j];
		}
	}
	return n;
}

vector<double> operator*(const Matrix& M, const vector<double>& vec) {
	if (M.columns != vec.size()) {
		throw "matrix*vector error";
	}
	vector<double> prod(M.rows, 0.0);
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			prod[i] += M.T[i][j] * vec[j];
		}
	}
	return prod;
}

Matrix operator*(const Matrix& M, double num) {
	Matrix new_M(M.rows, M.columns);
	for (int i = 0; i < M.rows; ++i) {
		for (int j = 0; j < M.columns; ++j) {
			new_M.T[i][j] = M.T[i][j] * num;
		}
	}
	return new_M;
}

class Assembler {
public:
	SparseMatrix* s_matrix;
	Grid* grid;

	Assembler(SparseMatrix* matrix) {
		s_matrix = matrix;
		grid = s_matrix->G;
	}


	int mu(int i) {
		return i % 3;
	}

	int nu(int j) {
		return (int)floor(j / 3.0);
	}

	double Gx(int i, int j, double h) {
		static Matrix M({{7, -8, 1},{-8, 16, -8},{1, -8, 7}});
		return M.T[i][j] / (3.0 * h);
	}

	double Gy(int i, int j, double h) {
		static Matrix M({ {7, -8, 1}, {-8, 16, -8}, {1, -8, 7} });
		return M.T[i][j] / (3.0 * h);
	}

	double Mx(int i, int j, double h) {
		static Matrix M({ {4,2,-1},{2,16,2},{-1,2,4} });
		return M.T[i][j] * (h / 30.0);
	}

	double My(int i, int j, double h) {
		static Matrix M({ {4,2,-1},{2,16,2},{-1,2,4} });
		return M.T[i][j] * (h / 30.0);
	}

	vector<double> get_f_vector() {
		vector<double> f;
		f.reserve(grid->nodes_number);
		for (int y = 0; y < grid->grid_y_size; ++y) {
			for (int x = 0; x < grid->grid_x_size; ++x) {
				f.push_back(get_f(grid->grid_x[x], grid->grid_y[y]));
			}
		}
		return f;
	}

	void put_matrix(const Matrix& matrix, const vector<int>& global_numbers) {
		for (int i = 0; i < global_numbers.size(); ++i) {
			for (int j = 0; j < i; ++j) {
				int start = min(0, s_matrix->ig[global_numbers[i]]);
				int end = max(s_matrix->ig[global_numbers[i] + 1] - s_matrix->ig[global_numbers[i]], (int)s_matrix->ig.size());
				auto iter = find(s_matrix->jg.begin() + start, s_matrix->jg.begin() + end, global_numbers[j]);
				s_matrix->gg[distance(s_matrix->jg.begin(), iter)] += matrix.T[i][j];
			}
			s_matrix->di[global_numbers[i]] += matrix.T[i][i];
		}
	}

	void put_vector(vector<double>& global_vector, const vector<double>& local_vector, const vector<int>& global_numbers) {
		for (int i = 0; i < local_vector.size(); ++i) {
			global_vector[global_numbers[i]] += local_vector[i];
		}
	}

	vector<double> Assemble(const vector<double>& global_f) {

		vector<double> global_b(global_f.size(), 0.0);

		for (auto& elem : grid->elements) {

			vector<double> f;
			f.reserve(9);
			for (int i : elem) {
				f.push_back(global_f[i]);
			}

			auto [first_x, first_y] = grid->number_to_coordinates[elem[0]];
			double hx = grid->grid_x[first_x + 2] - grid->grid_x[first_x];
			double hy = grid->grid_y[first_y + 2] - grid->grid_y[first_y];
			double lambda = get_lambda(grid->grid_x[first_x + 1], grid->grid_y[first_y + 1]);
			double gamma = get_gamma(grid->grid_x[first_x + 1], grid->grid_y[first_y + 1]);

			Matrix G(9, 9);
			Matrix M(9, 9);
			Matrix C(9, 9);
			
			for (int i = 0; i < 9; ++i) {
				for (int j = 0; j <= i; ++j) {
					double g = Gx(mu(i), mu(j), hx) * My(nu(i), nu(j), hy) + Mx(mu(i), mu(j), hx) * Gy(nu(i), nu(j), hy);
					g *= lambda;
					double c = Mx(mu(i), mu(j), hx) * My(nu(i), nu(j), hy);
					double m = c * gamma;
					G.T[i][j] = g;
					M.T[i][j] = m;
					C.T[i][j] = c;
					if (i != j) {
						G.T[j][i] = g;
						M.T[j][i] = m;
						C.T[j][i] = c;
					}
				}
			}

			auto A = M + G;
			auto b = C * f;

			put_matrix(A, elem);
			put_vector(global_b, b, elem);

		}

		return global_b;
	}

	void put_condition_three(const vector<int>& global_numbers, double betta, vector<double>& global_b) {

		auto [start_x, start_y] = grid->number_to_coordinates[global_numbers[0]];
		auto [end_x, end_y] = grid->number_to_coordinates[global_numbers[2]];

		double h = sqrt(pow(grid->grid_x[end_x] - grid->grid_x[start_x], 2) + pow(grid->grid_y[end_y] - grid->grid_y[start_y], 2));

		Matrix A({ {4, 2, -1},{2,16,2},{-1,2,4} });
		A = A * (betta * h / 30.0);
		vector<double> vec(3, 0.0);
		for (int i = 0; i < 3; ++i) {
			vec[i] = global_b[global_b[i]];
		}
		auto b = A * vec;
		put_matrix(A, global_numbers);
		put_vector(global_b, b, global_numbers);
	}

	void put_condition_two(const vector<int>& global_numbers, vector<double>& global_b) {
		auto [start_x, start_y] = grid->number_to_coordinates[global_numbers[0]];
		auto [end_x, end_y] = grid->number_to_coordinates[global_numbers[2]];
		double h = sqrt(pow(grid->grid_x[end_x] - grid->grid_x[start_x], 2) + pow(grid->grid_y[end_y] - grid->grid_y[start_y], 2));
		Matrix A({ {4, 2, -1},{2,16,2},{-1,2,4} });
		A = A * (h / 30.0);
		vector<double> vec(3, 0.0);
		for (int i = 0; i < 3; ++i) {
			vec[i] = 1;
		}
		auto b = A * vec;
		put_matrix(A, global_numbers);
		put_vector(global_b, b, global_numbers);
	}

	void put_condition_one(const vector<int>& global_numbers, vector<double>& global_b) {
		double num = 1.0e+20;
		for (int i = 0; i < global_numbers.size(); ++i) {
			s_matrix->di[global_numbers[i]] = num;
			auto [x, y] = grid->number_to_coordinates[global_numbers[i]];
			global_b[global_numbers[i]] = num * get_u(grid->grid_x[x], grid->grid_y[y]);
		}
	}

	//void put_condition_one(const vector<int>& global_numbers, vector<double>& global_b) {
	//	for (int i = 0; i < global_numbers.size(); ++i) {
	//		//cout << "-----\n";
	//		int number = global_numbers[i];
	//		auto [x, y] = grid->number_to_coordinates[number];
	//		double u = get_u(grid->grid_x[x], grid->grid_y[y]);
	//		global_b[number] = u;
	//		s_matrix->di[number] = 1.0;
	//
	//		for (int j = s_matrix->ig[number]; j < s_matrix->ig[number + 1]; ++j) {
	//			global_b[s_matrix->jg[j]] -= s_matrix->gg[j] * u;
	//			s_matrix->gg[j] = 0.0;
	//		}
	//
	//		for (int j = number + 1; j < s_matrix->G->nodes_number; ++j) {
	//			int start = min(0,s_matrix->ig[j]);
	//			int end = max(s_matrix->ig[j + 1] - s_matrix->ig[j], (int)s_matrix->ig.size());
	//			auto iter = find(s_matrix->jg.begin() + start, s_matrix->jg.begin() + end, number);
	//			if (iter == s_matrix->jg.end()) {
	//				continue;
	//			}
	//			int index = iter - s_matrix->jg.begin();
	//			global_b[j] -= s_matrix->gg[index] * u;
	//			s_matrix->gg[index] = 0.0;
	//			//cout << "!" << j << " \n";
	//		}
	//	}
	//}

	vector<double> get_u_vector() {
		vector<double> u;
		u.reserve(grid->nodes_number);
		for (int y = 0; y < grid->grid_y_size; ++y) {
			for (int x = 0; x < grid->grid_x_size; ++x) {
				u.push_back(get_u(grid->grid_x[x], grid->grid_y[y]));
			}
		}
		return u;
	}
};

