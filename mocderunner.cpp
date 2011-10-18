/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <string>
#include <stdio.h>

#include "mocde.h"
#include "benchmark.h"
#include "util.h"

using namespace std;

void showUsage(char *app) {
	printf("Usage:\n");
	printf("    %s function n\n", app);
}

int main(int argc, char **argv) {
	double F, CR, randomSeed;
	int maxEvaluations, populationSize, M;
	if (argc < 3) {
		showUsage(argv[0]);
		return EXIT_FAILURE;
	}
	
	int N = atoi(argv[2]);
	double bounds[N][2];
	char *instanceName = argv[1];
	benchmark::setup(instanceName, N, &M, bounds);
	
	bool silent = argc > 3 && string(argv[3]) == "--silent";
	if (!silent)
		printf("Population size: ");
	if (!scanf("%d", &populationSize)) {
		printf("Population size not given.");
		return EXIT_FAILURE;
	}
	if (!silent)
		printf("Maximum number of evaluations: ");
	if (!scanf("%d", &maxEvaluations)) {
		printf("Population size not given.");
		return EXIT_FAILURE;
	}
	if (!silent)
		printf("Differential variation (0..2): ");
	if (!scanf("%lf", &F)) {
		printf("Differential variation not given.");
		return EXIT_FAILURE;
	}
	if (!silent)
		printf("Crossover probability (0..1): ");
	if (!scanf("%lf", &CR)) {
		printf("Crossover probability not given.");
		return EXIT_FAILURE;
	}
	if (!silent)
		printf("Random seed: ");
	if (!scanf("%lf", &randomSeed)) {
		printf("Random seed not given.");
		return EXIT_FAILURE;
	}
	if (!silent)
		printf("Solving %s for %d parameters and %d objective functions.\n", instanceName, N, M);
	double **xs = util::createMatrix(populationSize, N);
	double **fxs = util::createMatrix(populationSize, M);
	
	MultiObjectiveCompactDifferentialEvolution de;
	int K = de.solve(xs, fxs, N, M, maxEvaluations, populationSize, CR, F, randomSeed, bounds, benchmark::evaluate);
	

	//printf("%s", util::toString(xs, K, N).c_str());
	printf("%s", util::toString(fxs, K, M).c_str());
	
	util::destroyMatrix(&xs, populationSize);
	util::destroyMatrix(&fxs, populationSize);

	return EXIT_SUCCESS;
}
