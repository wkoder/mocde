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
#include "util.h"
#include "moead/cec09.h"
#include "moead/global.h"

using namespace benchmark;

int nobj = 0;
int nreal = 0;
long rnd_uni_init = 90;
int seed = 177;
vector <double> idealpoint;
double scale[100];
double lowBound = 0;
double uppBound = 0;
int etax = 20;
int etam = 20;
double real;
double realm;
double realb = 0.9;

int evaluations = 0;
void (*function)(double *x, double *fx) = NULL;
double **xbounds;

bool OUT = false;

double **benchmark::getBounds() {
	return xbounds;
}

inline void setupBounds(double **bounds, int pos, double low, double up) {
	bounds[pos][0] = low;
	bounds[pos][1] = up;
}

void fsetup(void (*func)(double *x, double *fx), int real, int obj) {
	function = func;
	nreal = real;
	nobj = obj;
}

void fsetup(void (*func)(double *x, double *fx), int real, int obj, double low, double up, double **bounds) {
	fsetup(func, real, obj);
	xbounds = bounds;
	
	for (int i = 0; i < nreal; i++)
		setupBounds(bounds, i, low, up);
}

void fsetup(void (*func)(double *x, double *fx), int real, int obj, double fromlow, double tolow, double fromup, double toup, double **bounds) {
	fsetup(func, real, obj);
	xbounds = bounds;
	
	cerr << "Setting wide bounds for the MOEA/D implementation.\n";
	bounds[0][0] = fromlow;
	bounds[0][1] = fromup;
	double dlow = (tolow - fromlow) / (nreal-1);
	double dup = (toup - fromup) / (nreal-1);
	for (int i = 1; i < nreal; i++)
		setupBounds(bounds, i, bounds[i-1][0] + dlow, bounds[i-1][1] + dup);
}

void benchmark::setup(char *functionName, int real, int *obj, double **bounds) {
	evaluations = 0;
	if (strcmp(functionName, "deb2") == 0)
		fsetup(deb2, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "deb3") == 0)
		fsetup(deb3, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "fonseca2") == 0)
		fsetup(fonseca2, real, 2, -4, 4, bounds);
	else if (strcmp(functionName, "kursawe") == 0)
		fsetup(kursawe, real, 2, -5, 5, bounds);
	else if (strcmp(functionName, "wfg1") == 0)
		fsetup(wfg1, real, 2, 0, 0, 2, 2*real, bounds);
	else if (strcmp(functionName, "wfg2") == 0)
		fsetup(wfg2, real, 2, 0, 0, 2, 2*real, bounds);
	else if (strcmp(functionName, "wfg6") == 0)
		fsetup(wfg6, real, 2, 0, 0, 2, 2*real, bounds);
	else if (strcmp(functionName, "dtlz1") == 0)
		fsetup(dtlz1, real, 3, 0, 1, bounds);
	else if (strcmp(functionName, "dtlz2") == 0)
		fsetup(dtlz2, real, 3, 0, 1, bounds);
	else if (strcmp(functionName, "r_dtlz2") == 0)
		fsetup(r_dtlz2, real, 3, 0, 1, bounds);
	else if (strcmp(functionName, "dtlz3") == 0)
		fsetup(dtlz3, real, 3, 0, 1, bounds);
	else if (strcmp(functionName, "dtlz5im") == 0)
		fsetup(dtlz5im, real, 3, 0, 100, bounds);
	else if (strcmp(functionName, "dtlz7") == 0)
		fsetup(dtlz7, real, 3, 0, 100, bounds);
	else if (strcmp(functionName, "zdt1") == 0)
		fsetup(zdt1, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "zdt2") == 0)
		fsetup(zdt2, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "zdt3") == 0)
		fsetup(zdt3, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "uf1") == 0) {
		fsetup(CEC09::UF1, real, 2, -1, 1, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf2") == 0) {
		fsetup(CEC09::UF2, real, 2, -1, 1, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf3") == 0)
		fsetup(CEC09::UF3, real, 2, 0, 1, bounds);
	else if (strcmp(functionName, "uf4") == 0) {
		fsetup(CEC09::UF4, real, 2, -2, 2, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf5") == 0) {
		fsetup(CEC09::UF5, real, 2, -1, 1, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf6") == 0) {
		fsetup(CEC09::UF6, real, 2, -1, 1, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf7") == 0) {
		fsetup(CEC09::UF7, real, 2, -1, 1, bounds);
		setupBounds(bounds, 0, 0, 1);
	} else if (strcmp(functionName, "uf8") == 0) {
		fsetup(CEC09::UF8, real, 3, -2, 2, bounds);
		setupBounds(bounds, 0, 0, 1);
		setupBounds(bounds, 1, 0, 1);
	} else if (strcmp(functionName, "uf9") == 0) {
		fsetup(CEC09::UF9, real, 3, -2, 2, bounds);
		setupBounds(bounds, 0, 0, 1);
		setupBounds(bounds, 1, 0, 1);
	} else if (strcmp(functionName, "uf10") == 0) {
		fsetup(CEC09::UF10, real, 3, -2, 2, bounds);
		setupBounds(bounds, 0, 0, 1);
		setupBounds(bounds, 1, 0, 1);
	} else {
		printf("Benchmark instance %s not found.\n", functionName);
		exit(EXIT_FAILURE);
	}
	
	*obj = nobj;
}

void benchmark::evaluate(double *x, double *fx) {
	evaluations++;

//	double x2[nreal];
//	for (int i = 0; i < nreal; i++) {
//		if (x[i] < 0 || x[i] > 1)
//			cerr << "SHITTTTTTTT! Value: " << x[i] << endl;
//		x2[i] = xbounds[i][0] + x[i]*(xbounds[i][1] - xbounds[i][0]);
//	}
//	function(x2, fx);
	
//	bool out01 = false;
	bool canOut01 = false;
	for (int i = 0; i < nreal; i++) {
		if (x[i] < xbounds[i][0] || x[i] > xbounds[i][1])
			cerr << "SHITTTTTTTT! Value: " << x[i] << " out of [" << xbounds[i][0] << "," << xbounds[i][1] << "] at evalutation #" << evaluations << endl;
		if (xbounds[i][0] < 0 || xbounds[i][1] > 1) {
			canOut01 = true;
			if (x[i] < 0 || x[i] > 1) {
				OUT = true;
//				out01 = true;
//				cerr << i << ": " << x[i] << endl;
			}
		}
	}
//	if (canOut01 && !out01)
	if (evaluations == 100000 && canOut01 && !OUT) {
		cerr << "WARNING: Generating values only in [0..1] when in evaluation #" << evaluations << endl;
	}
//	if (evaluations % 1000 == 0)
//		cerr << evaluations << " " << canOut01 << " " << OUT << endl;
	
	function(x, fx);
}

int benchmark::getEvaluations() {
	return evaluations;
}
