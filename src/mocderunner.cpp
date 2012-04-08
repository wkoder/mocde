/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include "config.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "benchmark.h"
#include "util.h"
#include "randomlib.h"

#include "mocde.h"
#include "moead/algorithm.h"
#include "paes/paes.h"
#include "nsga2/global.h"
#include "mymoead.h"

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

int solve(double **xs, double **fxs, int n, int m, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx)) {
	randomize(randomSeed);
	initrandom(randomSeed * (1 << 30));
	
#ifdef MOCDE_IMPL
	MultiObjectiveCompactDifferentialEvolution de;
	return de.solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds, benchmark::evaluate);
#endif
#ifdef PAES_IMPL
	return paes(xs, fxs, bounds, function, 4, nreal, nobj, nreal*30, 2, populationSize, maxEvaluations, 0, 0.01, (int)(randomSeed*(1<<30)));
#endif
#ifdef NSGA2_IMPL
	return nsga2(xs, fxs, populationSize, maxEvaluations/populationSize, 0.9, 0.01, nreal, nobj, bounds, randomSeed);
#endif
#ifdef MOEAD_IMPL
	seed = (int)(randomSeed * (1 << 30));
	rnd_uni_init = 90.0;
	lowBound = 0;
	uppBound = 1;
	double **L = util::createMatrix(populationSize, nobj);
	int nicheSize = nobj == 2 ? 100 : 150;
	int updateLimit = nicheSize / 10;
		
	char filename[1024];
	sprintf(filename, "resources/W%dD.dat", nobj);
	ifstream file(filename);
	if (!file.is_open()) {
		cout << "File " << filename << " not found.\n";
		exit(-1);
	}
	for (int i = 0; i < populationSize; i++)
		for (int j = 0; j < nobj; j++)
			file >> L[i][j];
	file.close();
	
	CMOEAD MOEAD;
	MOEAD.load_parameter(populationSize, maxEvaluations/populationSize, nicheSize, updateLimit, F, maxEvaluations);
	return MOEAD.exec_emo(xs, fxs, L);
#endif
#ifdef MY_MOEAD_IMPL
	MyMOEAD moead;
	return moead.solve(xs, fxs, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds);
#endif
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
	if (populationSize <= 0) {
		cout << "Population size must be positive" << endl;
		exit(1);
	}
	if (!silent)
		cout << "Maximum number of evaluations: ";
	cin >> maxEvaluations;
	if (maxEvaluations <= 0) {
		cout << "Maximum number of evaluations must be positive" << endl;
		exit(1);
	}
	if (!silent)
		cout << "Differential variation [0..2]: ";
	cin >> F;
	if (F < 0 || F > 2) {
		cout << "Differential variation must be in [0..2]" << endl;
		exit(1);
	}
	if (!silent)
		cout << "Crossover probability [0..1): ";
	cin >> CR;
	if (CR < 0 || CR >= 1) {
		cout << "Crossover probability must be in [0..1)" << endl;
		exit(1);
	}
	if (!silent)
		cout << "Random seed [0..1): ";
	cin >> randomSeed;
	if (randomSeed < 0 || randomSeed >= 1) {
		cout << "Random seed must be in [0..1)" << endl;
		exit(1);
	}
	if (!silent)
		cout << "Solving " << instanceName << " for " << nreal << " parameters and " << nobj << " objective functions.\n";
	double **xs = util::createMatrix(populationSize, nreal);
	double **fxs = util::createMatrix(populationSize, nobj);
	
	int K = solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds, benchmark::evaluate);

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
