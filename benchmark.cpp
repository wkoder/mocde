/*
 * benchmark.c
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "benchmark.h"

#define INF 1e20

int instance = -1;
int n = 0;
int evaluations;
double **A;
double **B;
double *O;

double u(double x, double a, double b, double c) {
	if (x > a)
		return b * pow((x - a), c);
	if (x < a)
		return b * pow((-x - a), c);
	return 0;
}

double y(double x) {
	return 1 + (x + 1)/4;
}

void multiplyScalar(double *x, double c, double *r) {
	for (int i = 0; i < n; i++)
		r[i] = c * x[i];
}

double multiplyVector(double *x, double *y) {
	double f = 0;
	for (int i = 0; i < n; i++)
		f += x[i] * y[i];
	
	return f;
}

void subtractVector(double *x, double *y, double *r) {
	for (int i = 0; i < n; i++)
		r[i] = x[i] - y[i];
}

void multiplyMatrix(double *x, double **M, double *x2) {
	for (int i = 0; i < n; i++) {
		x2[i] = 0;
		for (int j = 0; j < n; j++)
			x2[i] += x[j] * M[j][i];
	}
}

// Sphere function
double sphere(double *x) {
	double f = 0;
	for (int i = 0; i < n; i++)
		f += x[i] * x[i];
	
	return f;
}

double rastrigins(double *z) {
	double f = 0;
	for (int i = 0; i < n; i++)
		f += z[i]*z[i] - 10*cos(2 * M_PI * z[i]) + 10;
	
	return f;
}

double weierstrass(double *z) {
	double a = 0.5;
	double b = 3.0;
	int kmax = 20;
	
	double f = 0;
	for (int i = 0; i < n; i++)
		for (int k = 0; k <= kmax; k++)
			f += pow(a, k) * cos(2 * M_PI * pow(b, k) * (z[i] + 0.5));
	
	return f;
}

double griewanks(double *x) {
	double f = 0;
	for (int i = 0; i < n; i++)
		f += x[i]*x[i] / 4000;
	
	return f;
}

double ackleys(double *x) {
	double f = 0;
	double a = 0;
	double b = 0;
	for (int i = 0; i < n; i++) {
		a += x[i] * x[i];
		b += cos(2 * M_PI * x[i]);
	}
	f = -20 * exp(-0.2 * sqrt(a / n)) - exp(b / n) + 20 + M_E;
	
	return f;
}

// Function composed by other 10 functions
double compose(double C, double *lambda, double *d, double *x, double ***M, double (**func)(double *x)) {
	double w[10];
	double fit[10];
	double z[n];
	double xl[n];
	double y[n];
	double yl[n];
	
	for (int i = 0; i < 10; i++)
		y[i] = 5;
	
	for (int i = 0; i < 10; i++) {
		double sum = sphere(x);
		w[i] = exp(-sum / (2*n*d[i]));
		
		multiplyScalar(x, 1/lambda[i], xl);
		multiplyMatrix(xl, M[i], z);
		fit[i] = func[i](z);
		
		multiplyScalar(y, 1/lambda[i], yl);
		multiplyMatrix(yl, M[i], z);
		double fmax = func[i](z);
		
		fit[i] = C * fit[i] / fmax;
	}
	
	double sw = 0;
	double maxw = 0;
	for (int i = 0; i < 10; i++) {
		sw += w[i];
		maxw = max(maxw, w[i]);
	}
	for (int i = 0; i < 10; i++) {
		if (w[i] != maxw)
			w[i] = w[i] * (1 - pow(maxw, 10));
		w[i] /= sw;
	}
	
	double f = 0;
	for (int i = 0; i < 10; i++)
		f += w[i] * fit[i];
	
	return f;
}

/**
 * References:
 * 	[42] A. K. Qin, V. L. Huang, and P. N. Suganthan, “Differential evolu-
		tion algorithm with strategy adaptation for global numerical optimiza-
		tion,” IEEE Trans. Evol. Comput., vol. 13, no. 2, pp. 398–417, Apr.
		2009.
	[43] P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y.-P. Chen, A. Auger,
		and S. Tiwari, “Problem definitions and evaluation criteria for the
		CEC 2005 special session on real-parameter optimization,” Nanyang
		Technol. Univ. KanGAL, Singapore, IIT Kanpur, Kanpur, India, Tech.
		Rep. 2 005 005, 2005.
	[44] J. Liang, P. Suganthan, and K. Deb, “Novel composition test functions
		for numerical global optimization,” in Proc. IEEE Symp. Swarm Intell.,
		2005, pp. 68–75.
	[45] J. Vesterstrøm and R. Thomsen, “A comparative study of differential
		evolution particle swarm optimization and evolutionary algorithms on
		numerical benchmark problems,” in Proc. IEEE Congr. Evol. Comput.,
		vol. 3. Jun. 2004, pp. 1980–1987.
 */

// Shifted sphere function: F1 from [43]
double f1(double *x) {
	double z[n];
	subtractVector(x, O, z);
	
	return sphere(z);
}

// Shifted Schwefel’s problem 1.2: F2 from [43]
double f2(double *x) {
	double f = 0;
	double z[n];
	subtractVector(x, O, z);
	for (int i = 0; i < n; i++)
		for (int j = 0; j <= i; j++)
			f += z[j] * z[j];
	
	return f;
}

// Rosenbrock’s function: f3 from [42]
double f3(double *x) {
	double f = 0;
	for (int i = 0; i < n-1; i++)
		f += 100 * pow(x[i]*x[i] - x[i+1], 2) + pow(x[i]-1, 2);
	
	return f;
}

// Shifted Ackley’s function: f5 from [42]
double f4(double *x) {
	return ackleys(x);
}

// Shifted rotated Ackley’s function: f6 from [42]
double f5(double *x) {
	double z[n];
	multiplyMatrix(x, A, z);
	
	return f4(z);
}

// Shifted Griewank’s function: f7 from [42]
double f6(double *x) {
	return griewanks(x);
}

// Shifted rotated Griewank’s function: f8 from [42]
double f7(double *x) {
	double z[n];
	multiplyMatrix(x, A, z);
	
	return f6(z);
}

// Shifted Rastrigin’s function: F9 from [43]
double f8(double *x) {
	double z[n];
	subtractVector(x, O, z);
	
	return rastrigins(z);
}

// Shifted rotated Rastrigin’s function: F10 from [43]
double f9(double *x) {
	double f = 0;
	double z[n];
	double tmp[n];
	subtractVector(x, O, tmp);
	multiplyMatrix(tmp, A, z);
	
	for (int i = 0; i < n; i++)
		f += z[i]*z[i] - 10*cos(2 * M_PI * z[i]) + 10;
	
	return f;
}

// Shifted noncontinuous Rastrigin’s function: f11 from [42]
double f10(double *x) {
	double f = 0;
	for (int i = 0; i < n; i++) {
		double y = fabs(x[i]) < 0.5 ? x[i] : round(2*x[i])/2;
		f += y*y - 10*cos(2 * M_PI * y) + 10;
	}
	
	return f;
}

// Schwefel’s function: f12 from [42]
double f11(double *x) {
	double f = 0;
	f += n * 418.9829;
	for (int i = 0; i < n; i++)
		f -= x[i] * sin(sqrt(fabs(x[i])));
	
	return f;
}

/**
 * Composition function 1: CF1 from [44]. 
 * The function f12 (CF1) is composed using ten sphere functions
 */
double f12(double *x) {
	double C = 2000;
	double lambda[10];
	double d[10];
	double **M[10];
	double ((*func[10])(double *x));
	
	for (int i = 0; i < 10; i++) {
		lambda[i] = 5 / 100.0;
		d[i] = 1;
		M[i] = A;
		func[i] = sphere;
	}
	
	return compose(C, lambda, d, x, M, func);
}

/**
 * Composition function 6: CF6 from [44]. 
 * The function f13 (CF6) is composed by using ten different benchmark functions, i.e., 
 * two rotated Rastrigin’s functions, two rotated Weierstrass functions, 
 * two rotated Griewank’s functions, two rotated Ackley’s functions, and two rotated Sphere functions.
 */
double f13(double *x) {
	double C = 2000;
	double lambda[10];
	double d[10];
	double **M[10];
	double ((*func[10])(double *x));
	
	for (int i = 0; i < 10; i++) {
		lambda[i] = 0.1 * (i+1);
		d[i] = 0.1 * (i+1);
		M[i] = A;
	}
	
	lambda[0] /= 5;
	lambda[1] /= 5;
	lambda[2] /= 0.5;
	lambda[3] /= 0.5;
	lambda[4] /= 100;
	lambda[5] /= 100;
	lambda[6] /= 32;
	lambda[7] /= 32;
	lambda[8] /= 100;
	lambda[9] /= 100;
	
	func[0] = func[1] = rastrigins;
	func[2] = func[3] = weierstrass;
	func[4] = func[5] = griewanks;
	func[6] = func[7] = ackleys;
	func[8] = func[9] = sphere;
	
	return compose(C, lambda, d, x, M, func);
}

// Schwefel problem 2.22: f2 from [45]
double f14(double *x) {
	double f = 0;
	double r = 1;
	for (int i = 0; i < n; i++) {
		f += fabs(x[i]);
		r *= fabs(x[i]);
	}
	f += r;
	
	return f;
}

// Schwefel problem 2.21: f4 from [45]
double f15(double *x) {
	double f = 0;
	for (int i = 0; i < n; i++)
		f = max(f, fabs(x[i]));
	
	return f;
}

// Generalized penalized function 1: f12 from [45]
double f16(double *x) {
	double f = 0;
	f += 10 * pow(sin(M_PI * y(x[0])), 2);
	for (int i = 0; i < n-1; i++)
		f += pow(y(x[i])-1, 2) * (1 + 10*pow(sin(M_PI * y(x[i+1])), 2));
	f += pow(y(x[n-1]) - 1, 2);
	f *= M_PI / n;
	for (int i = 0; i < n; i++)
		f += u(x[i], 10, 100, 4);
	
	return f;
}

// Generalized penalized function 2: f13 from [45]
double f17(double *x) {
	double f = 0;
	f += pow(sin(3 * M_PI * x[0]), 2);
	for (int i = 0; i < n-1; i++)
		f += pow(x[i]-1, 2) * (1 + pow(sin(3 * M_PI * x[i+1]), 2));
	f += (x[n-1]-1) * (1 + pow(sin(2 * M_PI * x[n-1]), 2));
	f *= 0.1;
	for (int i = 0; i < n; i++)
		f += u(x[i], 5, 100, 4);
	
	return f;
}

// Schwefel’s problem 2.6 with Global Optimum on Bounds: F5 from [43]
double f18(double *x) {
	double f = -INF;
	for (int i = 0; i < n; i++) {
		double r = multiplyVector(A[i], x) - multiplyVector(A[i], O);
		f = max(f, fabs(r));
	}
	
	return f;
}

// Shifted rotated Weierstrass function: F11 from [43]
double f19(double *x) {
	double f = 0;
	double a = 0.5;
	double b = 3.0;
	double z[n];
	double tmp[n];
	int kmax = 20;
	
	subtractVector(x, O, tmp);
	multiplyMatrix(tmp, A, z);
	f = weierstrass(z);
	for (int k = 0; k <= kmax; k++)
		f -= n * pow(a, k) * cos(2 * M_PI * pow(b, k) * 0.5);
	
	return f;
}

// Schwefel’s problem 2.13: F12 from [43]
double f20(double *x) {
	double f = 0;
	for (int i = 0; i < n; i++) {
		double a = 0;
		double b = 0;
		for (int j = 0; j < n; j++) {
			a += A[i][j]*sin(O[j]) + B[i][j]*cos(O[j]);
			b += A[i][j]*sin(x[j]) + B[i][j]*cos(x[j]);
		}
		
		f += (a-b) * (a-b);
	}
	
	return f;
}

void readInput(bool readO, bool readA, bool readB) {
	char fileName[1024];
	snprintf(fileName, sizeof(fileName), "input_data/f%d_%d.in", instance, n);
	
	FILE *fpt;
	fpt = fopen(fileName, "r");
	if (fpt == NULL) {
		printf("ERROR: Cannot open file %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	
	if (readO) {
		O = new double[n];
		for (int i = 0; i < n; i++)
			if (!fscanf(fpt, "%lf", &O[i])) {
				printf("ERROR: Cannot read V data from %s\n", fileName);
				exit(EXIT_FAILURE);
			}
	}
				
	if (readA) {
		A = new double*[n];
		for (int i = 0; i < n; i++)
			A[i] = new double[n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (!fscanf(fpt, "%lf", &A[i][j])) {
					printf("ERROR: Cannot read A data from %s\n", fileName);
					exit(EXIT_FAILURE);
				}
	}
	
	if (readB) {
		B = new double*[n];
		for (int i = 0; i < n; i++)
			B[i] = new double[n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (!fscanf(fpt, "%lf", &B[i][j])) {
					printf("ERROR: Cannot read B data from %s\n", fileName);
					exit(EXIT_FAILURE);
				}
	}
	
	fclose(fpt);
}

void setupBounds(double low, double up, double (*bounds)[2]) {
	for (int i = 0; i < n; i++) {
		bounds[i][0] = low;
		bounds[i][1] = up;
	}
}

void benchmarkSetupInstance(int inst, int nn, double (*bounds)[2]) {
	evaluations = 0;
	instance = inst;
	n = nn;
	A = NULL;
	B = NULL;
	O = NULL;
	
	switch (instance) {
		case 1: setupBounds(-100, 100, bounds); readInput(true, false, false); return;
		case 2: setupBounds(-100, 100, bounds); readInput(true, false, false); return;
		case 3: setupBounds(-100, 100, bounds); return;
		case 4: setupBounds(-32, 32, bounds); return;
		case 5: setupBounds(0, 600, bounds); readInput(false, true, false); return;
		case 6: setupBounds(0, 600, bounds); return;
		case 7: setupBounds(-5, 5, bounds); readInput(false, true, false); return;
		case 8: setupBounds(-5, 5, bounds); readInput(true, false, false); return;
		case 9: setupBounds(-5, 5, bounds); readInput(true, true, false); return;
		case 10: setupBounds(-500, 500, bounds); return;
		case 11: setupBounds(-500, 500, bounds); return;
		case 12: setupBounds(-5, 5, bounds); readInput(false, true, false); return;
		case 13: setupBounds(-5, 5, bounds); readInput(false, true, false); return;
		case 14: setupBounds(-10, 10, bounds); return;
		case 15: setupBounds(-100, 100, bounds); return;
		case 16: setupBounds(-50, 50, bounds); return;
		case 17: setupBounds(-50, 50, bounds); return;
		case 18: setupBounds(-100, 100, bounds); readInput(true, true, false); return;
		case 19: setupBounds(-0.5, 0.5, bounds); readInput(true, true, false); return;
		case 20: setupBounds(-M_PI, M_PI, bounds); readInput(true, true, true); return;
		default:
			printf("Benchmark instance %d not found.\n", instance);
			exit(EXIT_FAILURE);
	}
}

double benchmarkEvaluation(double *x) {
	evaluations++;
	switch (instance) {
		case 1:
			return f1(x);
		case 2:
			return f2(x);
		case 3:
			return f3(x);
		case 4:
			return f4(x);
		case 5:
			return f5(x);
		case 6:
			return f6(x);
		case 7:
			return f7(x);
		case 8:
			return f8(x);
		case 9:
			return f9(x);
		case 10:
			return f10(x);
		case 11:
			return f11(x);
		case 12:
			return f12(x);
		case 13:
			return f13(x);
		case 14:
			return f14(x);
		case 15:
			return f15(x);
		case 16:
			return f16(x);
		case 17:
			return f17(x);
		case 18:
			return f18(x);
		case 19:
			return f19(x);
		case 20:
			return f20(x);
		default:
			printf("Benchmark instance %d not found.\n", instance);
			exit(EXIT_FAILURE);
	}
}

int benchmarkGetEvaluations() {
	return evaluations;
}
