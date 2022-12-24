#pragma once

#include "grid.hpp"

#include <algorithm>
#include <cmath>

using namespace std;

class SparseMatrix {
public:
	Grid* G;
	int n;
	vector<double> di;
	vector<double> gg;
	vector<int> jg;
	vector<int> ig;

	SparseMatrix() {
		G = nullptr;
		n = 0;
		di.reserve(1);
		gg.reserve(1);
		ig.reserve(1);
		jg.reserve(1);
	}

	SparseMatrix(Grid* grid) {
		G = grid;
		n = G->nodes_number;
		di.resize(n, 0.0);
		jg.reserve(G->connections);
		ig.reserve(G->nodes_number + 1);
		ig.push_back(0);
		int counter = 0;
		for (auto adj : grid->adjacency_list) {
			for (int x : adj) {
				jg.push_back(x);
			}
			counter += adj.size();
			ig.push_back(counter);
		}
		gg.resize(G->connections, 0.0);
	}

	SparseMatrix Decomposition() {
		SparseMatrix matrix;
		matrix.n = n;
		matrix.di.assign(di.begin(), di.end());
		matrix.gg.assign(gg.begin(), gg.end());
		matrix.ig.assign(ig.begin(), ig.end());
		matrix.jg.assign(jg.begin(), jg.end());
		matrix.G = G;

		for (int i = 0; i < matrix.n; ++i) {
			double sum_d = 0.0;
			for (int j = 0; j < matrix.ig[i + 1]; ++j) {
				double sum_i_prev = 0.0;
				for (int k = matrix.ig[i]; k < j; ++k) {
					int i_prev = i - matrix.jg[j];
					int start = min(0,matrix.ig[i - i_prev]);
					int end = max(matrix.ig[i - i_prev + 1] - matrix.ig[i - i_prev], (int)matrix.jg.size());
					int value = matrix.jg[k];
					auto iter = find(matrix.jg.begin() + start, matrix.jg.begin() + end, value);
					if (iter != matrix.jg.end()) {
						sum_i_prev += matrix.gg[k] * matrix.gg[distance(matrix.jg.begin(),iter)];
					}
				}
				matrix.gg[j] = (matrix.gg[j] - sum_i_prev) / matrix.di[matrix.jg[j]];
				sum_d += matrix.gg[j] * matrix.gg[j];
			}
			matrix.di[i] = sqrt(matrix.di[i] - sum_d);
		}

		return matrix;
	}
};

double operator*(const vector<double>& v1, const vector<double>& v2) {
	if (v1.size() != v2.size()) throw "scalar_product error: incompatible vectors";
	double sum = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		sum += v1[i] * v2[i];
	}
	return sum;
}

vector<double> operator+(const vector<double>& v1, const vector<double>& v2) {
	if (v1.size() != v2.size()) throw "vector+ error: incompatible vectors";
	vector<double> sum(v1.size(), 0.0);
	for (int i = 0; i < v1.size(); ++i) {
		sum[i] = v1[i] + v2[i];
	}
	return sum;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2) {
	if (v1.size() != v2.size()) throw "vector- error: incompatible vectors";
	vector<double> su(v1.size(), 0.0);
	for (int i = 0; i < v1.size(); ++i) {
		su[i] = v1[i] - v2[i];
	}
	return su;
}

vector<double> operator* (const vector<double>& v1, double num) {
	vector<double> prod(v1.size(), 0.0);
	for (int i = 0; i < v1.size(); ++i) {
		prod[i] = v1[i] * num;
	}
	return prod;
}

double norm(const vector<double>& v1) {
	double sum = 0.0;
	for (int i = 0; i < v1.size(); ++i) {
		sum += v1[i] * v1[i];
	}
	return sqrt(sum);
}

vector<double> operator*(const SparseMatrix& M, const vector<double>& input) {
	vector<double> output(input.size(), 0.0);
	for (int i = 0; i < M.n; ++i) {
		output[i] = M.di[i] * input[i];
		for (int j = M.ig[i]; j < M.ig[i + 1]; ++j) {
			output[i] += M.gg[j] * input[M.jg[j]];
			output[M.jg[j]] += M.gg[j] * input[i];
		}
	}
	return output;
}

vector<double> solve_SLAE(const SparseMatrix& M, const vector<double>& b) {
	vector<double> x(M.n, 0.0);

	for (int i = 0; i < M.n; ++i) {
		double sum = 0.0;
		for (int j = M.ig[i]; j < M.ig[i + 1]; ++j) {
			sum += M.gg[j] * x[M.jg[j]];
		}
		x[i] = (b[i] - sum) / M.di[i];
	}

	for (int i = M.n - 1; i >= 0; i--) {
		x[i] /= M.di[i];
		for (int j = M.ig[i + 1] - 1; j >= M.ig[i]; j--) {
			x[M.jg[j]] -= M.gg[j] * x[i];
		}
	}

	return x;
}

vector<double> MSG_solve(SparseMatrix& A, const vector<double>& x0, const vector<double>& b, double eps, int max_iter) {

	auto choletsky_matrix = A.Decomposition();

	auto r = b - A * x0;
	auto z = solve_SLAE(choletsky_matrix, r);
	vector<double> x(x0.begin(), x0.end());
	double b_norm = norm(b);
	double residual = norm(r) / b_norm;
	
	for (int i = 1; i <= max_iter && residual > eps; ++i) {
		auto Mr = solve_SLAE(choletsky_matrix, r);
		double scalar_Mr_R = Mr * r;
		auto A_z = A * z;
		double alpha_k = scalar_Mr_R / (A_z * z);
		auto x_next = x + z * alpha_k;
		auto r_next = r - A_z * alpha_k;
		auto Mr_next = solve_SLAE(choletsky_matrix, r_next);
		double beta_k = (Mr_next * r_next) / scalar_Mr_R;
		auto z_next = solve_SLAE(choletsky_matrix, r_next) + z * beta_k;
		residual = norm(r_next) / b_norm;
		x = x_next;
		r = r_next;
		z = z_next;
	}

	return x;
}