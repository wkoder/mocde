/*
 * MultiObjectiveCompactDifferentialEvolution.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#ifndef MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_
#define MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_

class MultiObjectiveCompactDifferentialEvolution {
public:
	MultiObjectiveCompactDifferentialEvolution();
	virtual ~MultiObjectiveCompactDifferentialEvolution();
	
	int solve(double **xs, double **fxs, int n, int m, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx));

private:
	double solve(double *xs, int n, int maxEvaluations, int populationSize, double CR,
					double F, double randomSeed, double (*bounds)[2], double (*function)(double *x));
};

#endif /* MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_ */
