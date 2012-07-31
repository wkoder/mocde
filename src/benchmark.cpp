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
#include "util.h"
#include "randomlib.h"
#include "moead/cec09.h"
#include "moead/global.h"
#include "problems/deb.h"
#include "problems/misc.h"
#include "problems/zdt.h"
#include "problems/dtlz.h"
#include "problems/wfg.h"

using namespace benchmark;
using namespace CEC09;


int nobj = 0;
int nreal = 0;

vector <double> idealpoint;
double scale[100];
double lowBound = 0;
double uppBound = 0;

long rnd_uni_init;
int seed;
int etax = 20;
int etam = 20;
double real;
double realm;
double realb = 0.9;

int evaluations = 0;
void (*function)(double *x, double *fx) = NULL;
double **bounds;
double *varDelta;
double *xDelta;

bool OUT = false;

void checkBounds(double *x) {
	bool canOut01 = false;
	for (int i = 0; i < nreal; i++) {
		if (x[i] < bounds[i][0] || x[i] > bounds[i][1])
			cerr << "SHITTTTTTTT! Value: " << x[i] << " out of [" << bounds[i][0] << "," << bounds[i][1] << "] at evalutation #" << evaluations << endl;
		if (bounds[i][0] < 0 || bounds[i][1] > 1) {
			canOut01 = true;
			if (x[i] < 0 || x[i] > 1) {
				OUT = true;
			}
		}
	}
	
	if (evaluations % 5000 == 0 && canOut01 && !OUT) {
		cerr << "WARNING: Generating values only in [0..1] when in evaluation #" << evaluations << endl;
	}
}

inline void setupBounds(int pos, double low, double up) {
	bounds[pos][0] = low;
	bounds[pos][1] = up;
}

void fsetup(void (*func)(double *x, double *fx), int obj) {
	function = func;
	nobj = obj;
}

void fsetup(void (*func)(double *x, double *fx), int obj, double low, double up) {
	fsetup(func, obj);
	
	for (int i = 0; i < nreal; i++)
		setupBounds(i, low, up);
}

void fsetup(void (*func)(double *x, double *fx), int obj, double fromlow, double tolow, double fromup, double toup) {
	fsetup(func, obj);
	
	cerr << "Setting wide bounds for the MOEA/D implementation.\n";
	bounds[0][0] = fromlow;
	bounds[0][1] = fromup;
	double dlow = (tolow - fromlow) / (nreal-1);
	double dup = (toup - fromup) / (nreal-1);
	for (int i = 1; i < nreal; i++)
		setupBounds(i, bounds[i-1][0] + dlow, bounds[i-1][1] + dup);
}

double **benchmark::getBounds() {
	return bounds;
}

void benchmark::setup(char *functionName, int _nreal, int *_nobj) {
	nreal = _nreal;
	evaluations = 0;
	bounds = util::createMatrix(nreal, 2);
	xDelta = new double[nreal];
	
	if (strcmp(functionName, "deb2") == 0)
		fsetup(deb2, 2, 0, 1);
	else if (strcmp(functionName, "deb3") == 0)
		fsetup(deb3, 2, 0, 1);
	else if (strcmp(functionName, "fonseca2") == 0)
		fsetup(fonseca2, 2, -4, 4);
	else if (strcmp(functionName, "kursawe") == 0)
		fsetup(kursawe, 2, -5, 5);
	else if (strcmp(functionName, "wfg1") == 0)
		fsetup(wfg1, 2, 0, 0, 2, 2*nreal);
	else if (strcmp(functionName, "wfg2") == 0)
		fsetup(wfg2, 2, 0, 0, 2, 2*nreal);
	else if (strcmp(functionName, "wfg6") == 0)
		fsetup(wfg6, 2, 0, 0, 2, 2*nreal);
	else if (strcmp(functionName, "dtlz1") == 0)
		fsetup(dtlz1, 3, 0, 1);
	else if (strcmp(functionName, "dtlz2") == 0)
		fsetup(dtlz2, 3, 0, 1);
	else if (strcmp(functionName, "dtlz3") == 0)
		fsetup(dtlz3, 3, 0, 1);
	else if (strcmp(functionName, "dtlz4") == 0)
		fsetup(dtlz4, 3, 0, 1);
	else if (strcmp(functionName, "dtlz5") == 0)
		fsetup(dtlz5, 3, 0, 1);
	else if (strcmp(functionName, "dtlz6") == 0)
		fsetup(dtlz6, 3, 0, 1);
	else if (strcmp(functionName, "dtlz7") == 0)
		fsetup(dtlz7, 3, 0, 1);
//	else if (strcmp(functionName, "r_dtlz2") == 0)
//		fsetup(r_dtlz2, 3, 0, 1);
//	else if (strcmp(functionName, "dtlz5im") == 0)
//		fsetup(dtlz5im, 3, 0, 100);
	else if (strcmp(functionName, "zdt1") == 0)
		fsetup(zdt1, 2, 0, 1);
	else if (strcmp(functionName, "zdt2") == 0)
		fsetup(zdt2, 2, 0, 1);
	else if (strcmp(functionName, "zdt3") == 0)
		fsetup(zdt3, 2, 0, 1);
	else if (strcmp(functionName, "zdt4") == 0) {
		fsetup(zdt4, 2, -5, 5);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "zdt6") == 0)
		fsetup(zdt6, 2, 0, 1);
	else if (strcmp(functionName, "uf1") == 0) {
		fsetup(UF1, 2, -1, 1);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf2") == 0) {
		fsetup(UF2, 2, -1, 1);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf3") == 0)
		fsetup(UF3, 2, 0, 1);
	else if (strcmp(functionName, "uf4") == 0) {
		fsetup(UF4, 2, -2, 2);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf5") == 0) {
		fsetup(UF5, 2, -1, 1);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf6") == 0) {
		fsetup(UF6, 2, -1, 1);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf7") == 0) {
		fsetup(UF7, 2, -1, 1);
		setupBounds(0, 0, 1);
	} else if (strcmp(functionName, "uf8") == 0) {
		fsetup(UF8, 3, -2, 2);
		setupBounds(0, 0, 1);
		setupBounds(1, 0, 1);
	} else if (strcmp(functionName, "uf9") == 0) {
		fsetup(UF9, 3, -2, 2);
		setupBounds(0, 0, 1);
		setupBounds(1, 0, 1);
	} else if (strcmp(functionName, "uf10") == 0) {
		fsetup(UF10, 3, -2, 2);
		setupBounds(0, 0, 1);
		setupBounds(1, 0, 1);
	} else {
		printf("Benchmark instance %s not found.\n", functionName);
		exit(EXIT_FAILURE);
	}
	
	// Random delta, this way the ideal point changes to avoid bias
	varDelta = new double[nreal];
	for (int i = 0; i < nreal; i++) {
		varDelta[i] = rndreal(bounds[i][0], bounds[i][1]);
		bounds[i][0] += varDelta[i];
		bounds[i][1] += varDelta[i];
	}
	
	*_nobj = nobj;
}

void benchmark::evaluate(double *x, double *fx) {
	for (int i = 0; i < nreal; i++)
		xDelta[i] = x[i] - varDelta[i];
	
	evaluations++;
	checkBounds(x);
	function(xDelta, fx);
}

int benchmark::getEvaluations() {
	return evaluations;
}

void benchmark::destroy() {
	util::destroyMatrix(&bounds, nreal);
	delete [] xDelta;
	delete [] varDelta;
}

double *benchmark::getVariableDelta() {
	return varDelta;
}
