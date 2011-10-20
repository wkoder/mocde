/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>

#include "mocde.h"
#include "benchmark.h"
#include "util.h"

using namespace std;

void showUsage(char *app) {
	cout << "Usage:\n";
	cout << "    " << app << " function n varfile objfile\n";
	cout << "Where:\n";
	cout << "    function   Function name\n";
	cout << "    n          Number of variables\n";
	cout << "    varfile    File to output the Pareto set\n";
	cout << "    objfile    File to output the Pareto front\n";
}

int main(int argc, char **argv) {
	double F, CR, randomSeed;
	int maxEvaluations, populationSize, M;
	if (argc < 5) {
		showUsage(argv[0]);
		return EXIT_FAILURE;
	}
	
	int N = atoi(argv[2]);
	double bounds[N][2];
	char *instanceName = argv[1];
	benchmark::setup(instanceName, N, &M, bounds);
	
	bool silent = argc > 3 && string(argv[3]) == "--silent";
	if (!silent)
		cout << "Population size: ";
	cin >> populationSize;
	cout << "Population size not given.";
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
		cout << "Solving " << instanceName << " for " << N << " parameters and " << M << " objective functions.\n";
	double **xs = util::createMatrix(populationSize, N);
	double **fxs = util::createMatrix(populationSize, M);
	
	MultiObjectiveCompactDifferentialEvolution de;
	int K = de.solve(xs, fxs, N, M, maxEvaluations, populationSize, CR, F, randomSeed, bounds, benchmark::evaluate);

	ofstream psfile(argv[3]);
	if (psfile.is_open()) {
		psfile << util::toString(xs, K, N);
		psfile.close();
	} else
		cout << "Cannot open file " << argv[4];

	ofstream pffile(argv[4]);
	if (pffile.is_open()) {
		pffile << util::toString(fxs, K, M);
		pffile.close();
	} else
		cout << "Cannot open file " << argv[3];
	
	util::destroyMatrix(&xs, populationSize);
	util::destroyMatrix(&fxs, populationSize);

	return EXIT_SUCCESS;
}