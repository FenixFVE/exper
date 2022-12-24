#pragma once

#include <cmath>

double get_u(double x, double y) {
	return pow(x, 4) + pow(y, 4);
}
double get_lambda(double x, double y) {
	return 1.0;
}
double get_gamma(double x, double y) {
	return 0.0;
}
double get_f(double x, double y) {
	return -12 * x * x - 12 * y * y;
}
