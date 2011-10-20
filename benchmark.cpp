/*
 * benchmark.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstring>

#include "benchmark.h"
#include "problemdef.h"

using namespace benchmark;

#define INF 1e20

int evaluations = 0;
void (*function)(double *x, double *fx) = NULL;

void fsetup(void (*func)(double *x, double *fx), int real, int obj, double low, double up, double (*bounds)[2]) {
	function = func;
	nreal = real;
	nobj = obj;
	for (int i = 0; i < nreal; i++) {
		bounds[i][0] = low;
		bounds[i][1] = up;
	}
}

void benchmark::setup(char *functionName, int real, int *obj, double (*bounds)[2]) {
	evaluations = 0;
	if (strcmp(functionName, "wfg1") == 0)
		fsetup(wfg1, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "wfg2") == 0)
		fsetup(wfg2, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "wfg6") == 0)
		fsetup(wfg6, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz1") == 0)
		fsetup(dtlz1, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz2") == 0)
		fsetup(dtlz2, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "r_dtlz2") == 0)
		fsetup(r_dtlz2, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz3") == 0)
		fsetup(dtlz3, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz5im") == 0)
		fsetup(dtlz5im, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz7") == 0)
		fsetup(dtlz7, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "zdt1") == 0)
		fsetup(zdt1, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "zdt2") == 0)
		fsetup(zdt2, real, 2, 0, 100, bounds);
	else if (strcmp(functionName, "zdt3") == 0)
		fsetup(zdt3, real, 2, 0, 100, bounds);
//	else if (strcmp(functionName, "deb1") == 0)
//		fsetup(deb1, real, 2, 0, 100, bounds);
//	else if (strcmp(functionName, "deb2") == 0)
//		fsetup(deb2, real, 2, 0, 100, bounds);
	else {
		printf("Benchmark instance %s not found.\n", functionName);
		exit(EXIT_FAILURE);
	}

	*obj = nobj;
}

void benchmark::evaluate(double *x, double *fx) {
	evaluations++;
	function(x, fx);
}

int benchmark::getEvaluations() {
	return evaluations;
}
