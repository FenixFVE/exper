
#include "all.hpp"
#include "t.hpp"

#include <set>
#include <vector>
#include <iostream>
#include <iomanip>


using namespace std;

int main() {
	Grid grid("grid.txt");
	SparseMatrix matrix(&grid);
	Assembler assembler(&matrix);
	auto u_star = assembler.get_u_vector();
	auto f = assembler.get_f_vector();
	auto b = assembler.Assemble(f);

	//assembler.put_condition_three({ 0, 1, 2 }, 1, b);
	//assembler.put_condition_two({ 6, 7, 8 }, b);
	assembler.put_condition_one({ 0, 1, 2, 3, 4, 5, 6, 7, 8  }, b);
	
	vector<double> x0(f.size(), 0.0);
	auto q = MSG_solve(matrix, x0, b, 10e-16, 1000);
	int str = 15;
	int counter = 0;
	double minus_norm = 0.0, norm = 0.0;
	printf("X |Y |\t\tU |\t    U* |     U* - U\n");
	for (int y = 0; y < grid.grid_y_size; ++y) {
		for (int x = 0; x < grid.grid_x_size; ++x) {
			cout << grid.grid_x[x] << " |" << grid.grid_y[y] << " | ";
			printf("%10f | ", q[counter]); 
			printf("%10f | ", u_star[counter]);
			printf("%10f\n", u_star[counter] - q[counter]);
			minus_norm += pow(u_star[counter] - q[counter], 2);
			norm += u_star[counter] * u_star[counter];
			counter++;
		}
	}
	minus_norm = sqrt(minus_norm); norm = sqrt(norm);
	cout << "|U*-U|/|U*| = " << minus_norm / norm;

	cout << "\nFINISH\n";
	return 0;
}
