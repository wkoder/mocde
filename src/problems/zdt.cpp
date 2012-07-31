/* Test problem definitions */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "zdt.h"
#include "../moead/global.h"

/*  Test problem ZDT1
 # of real variables = 30
 # of bin variables = 0
 # of objectives = 2
 # of constraints = 0
 */
void zdt1(double *xreal, double *obj) {
	double f1, f2, g, h;
	int i;
	f1 = xreal[0];
	g = 0.0;
	for (i = 1; i < nreal; i++) {
		g += xreal[i];
	}
	g = 9.0 * g / 29.0;
	g += 1.0;
	h = 1.0 - sqrt(f1 / g);
	f2 = g * h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

/*  Test problem ZDT2
 # of real variables = 30
 # of bin variables = 0
 # of objectives = 2
 # of constraints = 0
 */
void zdt2(double *xreal, double *obj) {
	double f1, f2, g, h;
	int i;
	f1 = xreal[0];
	g = 0.0;
	for (i = 1; i < nreal; i++) {
		g += xreal[i];
	}
	g = 9.0 * g / 29.0;
	g += 1.0;
	h = 1.0 - pow((f1 / g), 2.0);
	f2 = g * h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

/*  Test problem ZDT3
 # of real variables = 30
 # of bin variables = 0
 # of objectives = 2
 # of constraints = 0
 */
void zdt3(double *xreal, double *obj) {
	double f1, f2, g, h;
	int i;
	f1 = xreal[0];
	g = 0.0;
	for (i = 1; i < nreal; i++) {
		g += xreal[i];
	}
	g = 9.0 * g / 29.0;
	g += 1.0;
	h = 1.0 - sqrt(f1 / g) - (f1 / g) * sin(10.0 * M_PI * f1);
	f2 = g * h;
	obj[0] = f1;
	obj[1] = f2;
	return;
}

/* ZDT4 - MOP */
void zdt4(double *xreal, double *obj) {
	double f1, f2, g, h, sum;
	int j;

	f1 = xreal[0];
	sum = 0.0;
	for (j = 1; j < nreal; j++) {
		sum += (pow(xreal[j], 2.0) - 10.0 * cos(4.0 * M_PI * xreal[j]));
	}
	g = 1.0 + 10.0 * (nreal - 1.0) + sum;
	h = 1.0 - sqrt(f1 / g);
	f2 = g * h;
	obj[0] = f1;
	obj[1] = f2;
}

/* ZDT6 - MOP */
void zdt6(double *xreal, double *obj) {
	double f1, f2, g, h, sum;
	int j;

	f1 = 1.0 - exp(-4.0 * xreal[0]) * pow(sin(6.0 * M_PI * xreal[0]), 6.0);
	sum = 0.0;
	for (j = 1; j < nreal; j++) {
		sum += xreal[j];
	}
	g = 1.0 + 9.0 * pow((sum / (nreal - 1.0)), 0.25);
	h = 1.0 - pow((f1 / g), 2.0);
	f2 = g * h;
	obj[0] = f1;
	obj[1] = f2;
	obj[1] = f2;
}
