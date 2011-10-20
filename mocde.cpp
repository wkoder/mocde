/*
 * MultiObjectiveCompactDifferentialEvolution.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <math.h>
#include <stdio.h>

#include "mocde.h"
#include "rand.h"
#include "randomlib.h"
#include "benchmark.h"
#include "util.h"

MultiObjectiveCompactDifferentialEvolution::MultiObjectiveCompactDifferentialEvolution() {

}

MultiObjectiveCompactDifferentialEvolution::~MultiObjectiveCompactDifferentialEvolution() {

}

void ensureBounds(double *x, double (*bounds)[2], int n) {
	for (int i = 0; i < n; i++) {
		x[i] = max(x[i], bounds[i][0]);
		x[i] = min(x[i], bounds[i][1]);
	}
}

void generateX(double *x, double *u, double *d, double (*bounds)[2], int n) {
	for (int i = 0; i < n; i++)
		x[i] = normal(u[i], d[i]);
	ensureBounds(x, bounds, n);
}

void printx(char *d, double *x, int n) {
	printf("%s = ", d);
	for (int i = 0; i < n; i++)
		printf("%.6f ", x[i]);
	printf("\n");
}

double MultiObjectiveCompactDifferentialEvolution::solve(double *xb, int n, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], double (*function)(double *x)) {
	double u[n];
	double d[n];
	double elite[n];
	double xr[n];
	double xs[n];
	double xt[n];
	double xoff[n];
	warmup_random(randomSeed);
	initrandom(randomSeed * (1 << 30));
	int startingEval = benchmark::getEvaluations();
	
	// PV initialization
	for (int i = 0; i < n; i++) {
		double r = bounds[i][1] - bounds[i][0];
		u[i] = bounds[i][0] + r/2;
		d[i] = r * 0.341;
	}
	
	generateX(elite, u, d, bounds, n);
	double felite = function(elite);
	while (benchmark::getEvaluations()-startingEval < maxEvaluations) {
//		if (benchmark::getEvaluations() % populationSize == 0) {
//			printf("Iteration #%d:\n", benchmark::getEvaluations() / populationSize);
//			printf("	best: %s\n", util::toString(elite, n).c_str());
//			printf("	f(best): %.6f\n", felite);
//		}
		
		// Mutation
		generateX(xr, u, d, bounds, n);
		generateX(xs, u, d, bounds, n);
		generateX(xt, u, d, bounds, n);
		
		for (int i = 0; i < n; i++)
			xoff[i] = xt[i] + F*(xr[i] - xs[i]);
		ensureBounds(xoff, bounds, n);
		
		// Crossover
		for (int i = 0; i < n; i++)
			if (!flip(CR))
				xoff[i] = elite[i];
		
		// Elite selection
		double fxoff = function(xoff);
		double *winner = elite;
		double *loser = xoff;
		if (fxoff < felite) {
			winner = xoff;
			loser = elite;
			felite = fxoff;
		}
		
		// PV update
		for (int i = 0; i < n; i++) {
			double u2 = u[i] + (winner[i] - loser[i]) / populationSize;
			double d2 = sqrt(fabs(d[i]*d[i] + u[i]*u[i] - u2*u2 + 
					(winner[i]*winner[i] - loser[i]*loser[i]) / populationSize));
			
			u2 = max(u2, bounds[i][0]);
			u2 = min(u2, bounds[i][1]);
			u[i] = u2;
			d[i] = d2;
		}
		
		if (winner == xoff)
			for (int i = 0; i < n; i++)
				elite[i] = xoff[i];
	}
	
//	printf("Iteration #%d:\n", benchmark::getEvaluations() / populationSize);
//	printf("	best: %s\n", util::toString(elite, n).c_str());
//	printf("	f(best): %.6f\n", felite);
	
	for (int i = 0; i < n; i++)
		xb[i] = elite[i];
	
	return felite;
}

namespace mo {
	void (*mof)(double *x, double *fx);
	int N, M;
	double *L, *fx;
	
	double chebyshev(double *x) {
		double sum = 0;
		mof(x, fx);
		for (int i = 0; i < N; i++)
			sum += fx[i] * L[i];
		return sum;
	}
}

int MultiObjectiveCompactDifferentialEvolution::solve(double **xb, double **fxb, int n, int m, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx)) {
	using namespace mo;
	N = n;
	M = m;
	mof = function;
	L = new double[M];
	fx = new double[M];

	double dx = 1.0 / (populationSize - 1);
	for (int i = 0; i < populationSize; i++) {
		L[0] = i * dx;
		L[1] = 1 - L[0];
		solve(xb[i], N, maxEvaluations/populationSize, populationSize, CR, F, randomSeed, bounds, chebyshev);
		function(xb[i], fxb[i]);
	}

	delete [] L;
	delete [] fx;
	return populationSize;
}
