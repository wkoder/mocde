/* Test problem definitions */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "misc.h"
#include "../moead/global.h"

void fonseca2(double *xreal, double *obj) {
	double s1, s2;
	s1 = s2 = 0.0;
	for (int i = 0; i < nreal; i++) {
		s1 += pow((xreal[i] - (1.0 / sqrt((double) nreal))), 2.0);
		s2 += pow((xreal[i] + (1.0 / sqrt((double) nreal))), 2.0);
	}

	obj[0] = 1.0 - exp(-s1);
	obj[1] = 1.0 - exp(-s2);
}

void kursawe(double *xreal, double *obj) {
	double res1, res2;
	res1 = -0.2 * sqrt((xreal[0] * xreal[0]) + (xreal[1] * xreal[1]));
	res2 = -0.2 * sqrt((xreal[1] * xreal[1]) + (xreal[2] * xreal[2]));
	obj[0] = -10.0 * (exp(res1) + exp(res2));
	obj[1] = 0.0;
	for (int i = 0; i < 3; i++)
		obj[1] += pow(fabs(xreal[i]), 0.8) + 5.0 * sin(pow(xreal[i], 3.0));
}
