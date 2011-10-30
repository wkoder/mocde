/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "mocde.h"
#include "benchmark.h"
#include "util.h"

using namespace std;

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

void showUsage(char *app) {
	cout << "Usage:\n";
	cout << "    " << app << " function n varfile objfile [--silent]\n";
	cout << "Where:\n";
	cout << "    function   Function name\n";
	cout << "    n          Number of variables\n";
	cout << "    varfile    File to output the Pareto set\n";
	cout << "    objfile    File to output the Pareto front\n";
	cout << "    --silent   Sink all output.\n";
}

int main(int argc, char **argv) {
	double F, CR, randomSeed;
	int maxEvaluations, populationSize, nobj;
	if (argc < 5) {
		showUsage(argv[0]);
		return EXIT_FAILURE;
	}
	
	int nreal = atoi(argv[2]);
	double bounds[nreal][2];
	char *instanceName = argv[1];
	benchmark::setup(instanceName, nreal, &nobj, bounds);
	
	bool silent = argc > 5 && string(argv[5]) == "--silent";
	if (!silent)
		cout << "Population size: ";
	cin >> populationSize;
	if (!silent)
		cout << "Maximum number of evaluations: ";
	cin >> maxEvaluations;
	if (!silent)
		cout << "Differential variation (0..2): ";
	cin >> F;
	if (!silent)
		cout << "Crossover probability (0..1): ";
	cin >> CR;
	if (!silent)
		cout << "Random seed: ";
	cin >> randomSeed;
	if (!silent)
		cout << "Solving " << instanceName << " for " << nreal << " parameters and " << nobj << " objective functions.\n";
	double **xs = util::createMatrix(populationSize, nreal);
	double **fxs = util::createMatrix(populationSize, nobj);
	
	MultiObjectiveCompactDifferentialEvolution de;
	int K = de.solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds, benchmark::evaluate);

	ofstream psfile(argv[3]);
	if (psfile.is_open()) {
		psfile << util::toString(xs, K, nreal);
		psfile.close();
	} else
		cout << "Cannot open file " << argv[3];

	ofstream pffile(argv[4]);
	if (pffile.is_open()) {
		pffile << util::toString(fxs, K, nobj);
		pffile.close();
	} else
		cout << "Cannot open file " << argv[4];
	
	cout << "Evaluations: " << benchmark::getEvaluations() << endl;
//	cout << "x* stats:\n";
//	printStats(xs, populationSize, nreal);
//	cout << "f(x*) stats:\n";
//	printStats(fxs, populationSize, nobj);
	
	util::destroyMatrix(&xs, populationSize);
	util::destroyMatrix(&fxs, populationSize);

	return EXIT_SUCCESS;
}
