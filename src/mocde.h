/*
 * MultiObjectiveCompactDifferentialEvolution.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include "config.h"

#ifdef MOCDE_IMPL

#ifndef MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_
#define MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_

#include <vector>
#include "individual.h"

class MultiObjectiveCompactDifferentialEvolution {
public:
	MultiObjectiveCompactDifferentialEvolution();
	virtual ~MultiObjectiveCompactDifferentialEvolution();
	
	int solve(double **xs, double **fxs, int n, int m, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double **bounds, void (*function)(double *x, double *fx));

	int nobj;
	void (*function)(double *x, double *fx);

private:
	double solve(double *xs, int n, int maxEvaluations, int populationSize, double CR,
					double F, double *startingMean, double *startingStdDev, double **bounds, double (*function)(double *x));
	
	bool addToArchive(std::vector<Individual *> &archive, Individual *ind, Individual *parent);
	
	int nreal;
	int populationSize;

	// New implementation
	void initMOEADArchive(std::vector<Individual *> &archive);
	bool addToMOEADArchive(std::vector<Individual *> &archive, Individual *ind, Individual *parent);
	double chebyshevScalarizing(double *fx, double *namda);
	double **L;
	double *ideal;
	void updateIdeal(double *fx);
};

#endif /* MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_ */
#endif
