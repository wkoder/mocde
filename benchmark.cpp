/*
 * benchmark.c
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstring>

#include "benchmark.h"

using namespace benchmark;

#define INF 1e20

int N = 0;
int M = 0;
int evaluations = 0;
void (*function)(double *x, double *fx) = NULL;

void a(double *f, double *fx) {
	for (int i = 0; i < N; i++)
		fx[i] = abs(f[i]);
}

void b(double *f, double *fx) {
	for (int i = 0; i < N; i++)
		fx[i] = f[i]*f[i];
}

void fsetup(void (*func)(double *x, double *fx), int n, int m, double low, double up, double (*bounds)[2]) {
	function = func;
	N = n;
	M = m;
	for (int i = 0; i < N; i++) {
		bounds[i][0] = low;
		bounds[i][1] = up;
	}
}

void benchmark::setup(char *functionName, int n, int *m, double (*bounds)[2]) {
	evaluations = 0;
	if (strcmp(functionName, "deb1") == 0)
		fsetup(a, n, n, -100, 100, bounds);
	else if (strcmp(functionName, "deb2") == 0)
		fsetup(b, n, n, -100, 100, bounds);
	else {
		printf("Benchmark instance %s not found.\n", functionName);
		exit(EXIT_FAILURE);
	}

	*m = M;
}

void benchmark::evaluate(double *x, double *fx) {
	evaluations++;
	function(x, fx);
}

int benchmark::getEvaluations() {
	return evaluations;
}
