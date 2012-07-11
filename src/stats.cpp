/*
 * stats.cpp
 *
 *  Created on: Jun 15, 2012
 *      Author: Moises Osorio
 */

#include "stats.h"

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "util.h"
#include "benchmark.h"

using namespace std;

namespace stats {

static int outputInterval = 1e9;
static string fileNamePrefix;
static int nextReport = 1e9;
static double **tempMatrix = NULL;

void configure(int outputInt, string filePrefix) {
	outputInterval = outputInt;
	nextReport = outputInterval;
	fileNamePrefix = filePrefix;
	tempMatrix = util::createMatrix(500, 200);
}

bool timeToReport() {
	return benchmark::getEvaluations() >= nextReport;
}

void report(vector<Individual *> population) {
	if (!timeToReport())
		return;
	
	string id = util::toString(benchmark::getEvaluations());
	util::getPopulationX(tempMatrix, population);
	util::writeMatrixFile(tempMatrix, population.size(), population[0]->nreal, fileNamePrefix + "_var_" + id + ".out");
	util::getPopulationFx(tempMatrix, population);
	util::writeMatrixFile(tempMatrix, population.size(), population[0]->nobj, fileNamePrefix + "_front_" + id + ".out");
	nextReport += outputInterval;
}


void printStats(double *x, int n, int width, int precision) {
	for (int i = 0; i < n; i++)
		cout << setw(width) << setiosflags(ios::fixed) << setprecision(precision) << x[i] << " ";
	cout << endl;
}

void printStats(double **x, int r, int c) {
	double avg[c];
	double std[c];
	double stdr[c];
	
	for (int j = 0; j < c; j++) {
		avg[j] = 0;
		for (int i = 0; i < r; i++)
			avg[j] += x[i][j];
		
		avg[j] /= r;
		std[j] = 0;
		for (int i = 0; i < r; i++)
			std[j] += (x[i][j] - avg[j]) * (x[i][j] - avg[j]);
		std[j] = sqrt(std[j] / r);
		stdr[j] = avg[j] < EPS ? 0 : std[j] * 100 / avg[j];
	}
	
	int precision = 3;
	int width = precision + 4;
	printStats(avg, c, width, precision);
	printStats(std, c, width, precision);
	printStats(stdr, c, width, precision);
	cout << endl;
}


} /* namespace stats */
