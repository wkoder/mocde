/*
 * mocde.cpp
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <string>
#include <stdio.h>

#include "compactdifferentialevolution.h"
#include "benchmark.h"
#include "util.h"

using namespace std;

int main(int argc, char **argv) {
	double F, CR, randomSeed;
	int maxEvaluations, populationSize;
	
	int n = atoi(argv[2]);
	double bounds[n][2];
	double xs[n];
	benchmarkSetupInstance(atoi(argv[1]), n, bounds);
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
	
	CompactDifferentialEvolution de;
	double Fxs = de.solve(xs, n, maxEvaluations, populationSize, CR, F, 
			randomSeed, bounds, benchmarkEvaluation);
	
	printf("x* = %s\n", toString(xs, n).c_str());
	printf("F(x*) = %e\n", Fxs);
	
	return EXIT_SUCCESS;
}
