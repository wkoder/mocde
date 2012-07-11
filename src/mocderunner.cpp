/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include "config.h"

#include <stdlib.h>
#include <string>

#include "benchmark.h"
#include "util.h"
#include "randomlib.h"
#include "stats.h"

#include "mocde.h"
#include "moead/algorithm.h"
#include "paes/paes.h"
#include "nsga2/global.h"
#include "mymoead.h"
#include "jmetal/PAES_main.h"

using namespace std;

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

#ifdef MOCDE_IMPL
int solve(double **xs, double **fxs, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
				double F, int maxSurvival, double randomSeed, double **bounds, void (*function)(double *x, double *fx)) {
#else
	int solve(double **xs, double **fxs, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
					double F, double randomSeed, double **bounds, void (*function)(double *x, double *fx)) {
#endif
	randomize(randomSeed);
	initrandom(randomSeed * (1 << 30));
	
#ifdef MOCDE_IMPL
	MultiObjectiveCompactDifferentialEvolution de;
	return de.solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, maxSurvival, randomSeed, bounds, benchmark::evaluate);
#endif
#ifdef PAES_IMPL
#ifdef JMETAL
	PAES_main paes;
	return paes.solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds);
#else
	return paes(xs, fxs, bounds, function, 4, nreal, nobj, nreal*30, 2, populationSize, maxEvaluations, 0, 0.01, (int)(randomSeed*(1<<30)));
#endif
#endif
#ifdef NSGA2_IMPL
	return nsga2(xs, fxs, populationSize, maxEvaluations/populationSize, 0.9, 0.01, nreal, nobj, bounds, randomSeed);
#endif
#ifdef MOEAD_IMPL
	seed = (int)(randomSeed * (1 << 30));
	rnd_uni_init = 90.0;
	int nicheSize = nobj == 2 ? 100 : 150;
	int updateLimit = nicheSize / 10;
	
	CMOEAD MOEAD;
	MOEAD.load_parameter(populationSize, maxEvaluations/populationSize, nicheSize, updateLimit, F, maxEvaluations, bounds);
	return MOEAD.exec_emo(xs, fxs);
#endif
#ifdef MY_MOEAD_IMPL
	MyMOEAD moead;
	return moead.solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, bounds);
#endif
}

int main(int argc, char **argv) {
	double F, CR, randomSeed;
	int maxEvaluations, populationSize, nobj, maxSurvival;
	if (argc < 4) {
		showUsage(argv[0]);
		return EXIT_FAILURE;
	}
	
	int nreal = atoi(argv[2]);
	char *instanceName = argv[1];
	benchmark::setup(instanceName, nreal, &nobj);
	
	bool silent = argc > 4 && string(argv[4]) == "--silent";
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
#ifdef MOCDE_IMPL
	if (!silent)
		cout << "Maximum survival: ";
	cin >> maxSurvival;
	if (maxSurvival <= 0) // Infinite
		maxSurvival = 1<<30;
#endif
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
	string filePrefix(argv[3]);
	
	stats::configure(5000, filePrefix);
#ifdef MOCDE_IMPL
	int K = solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, maxSurvival, randomSeed, benchmark::getBounds(), benchmark::evaluate);
#else
	int K = solve(xs, fxs, nreal, nobj, maxEvaluations, populationSize, CR, F, randomSeed, benchmark::getBounds(), benchmark::evaluate);
#endif
	double *varDelta = benchmark::getVariableDelta();
	for (int i = 0; i < populationSize; i++)
		for (int j = 0; j < nreal; j++)
			xs[i][j] -= varDelta[j];
	
	util::writeMatrixFile(xs, K, nreal, filePrefix + "_var.out");
	util::writeMatrixFile(fxs, K, nobj, filePrefix + "_front.out");
	
	if (maxEvaluations != benchmark::getEvaluations()) {
		int extra = benchmark::getEvaluations() - maxEvaluations;
		cerr << "Did " << benchmark::getEvaluations() << " evaluations instead of " << maxEvaluations 
				<< " (" << extra << ")" << endl;
		if (extra > 0)
			return 1001;
	}
//	cout << "x* stats:\n";
//	printStats(xs, populationSize, nreal);
//	cout << "f(x*) stats:\n";
//	printStats(fxs, populationSize, nobj);
	
	util::destroyMatrix(&xs, populationSize);
	util::destroyMatrix(&fxs, populationSize);
	benchmark::destroy();

	return EXIT_SUCCESS;
}
