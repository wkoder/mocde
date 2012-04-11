/* Test problem definitions */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "deb.h"
#include "../moead/global.h"

void deb2(double *xreal, double *obj) {
	double g;
	double h;
	double f1;

	f1 = xreal[0];
	g = 1 + 10 * xreal[1];
	h = 1 - pow((f1 / g), 2) - (f1 / g) * sin(12 * M_PI * f1);

	obj[0] = f1;
	obj[1] = g * h;
}

void deb3(double *xreal, double *obj) {
	double g;
	double h;
	double f1;
	double alfa;
	double beta;
	double seno;
	double pi = 3.1415926535;
	double argumento;
	double q;

	alfa = 10;
	q = 10;
	beta = 1;

	argumento = q * pi * xreal[0];
	seno = sin(argumento);
	f1 = 1 - exp(-4 * xreal[0]) * pow(seno, 4);
	g = 1 + xreal[1] * xreal[1];
	if (f1 <= beta * g) {
		h = 1 - pow((f1 / (beta * g)), alfa);
	} else {
		h = 0;
	}

	obj[0] = f1;
	obj[1] = g * h;
}
