/*
 * CompactDifferentialEvolution.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#ifndef COMPACTDIFFERENTIALEVOLUTION_H_
#define COMPACTDIFFERENTIALEVOLUTION_H_

class CompactDifferentialEvolution {
public:
	CompactDifferentialEvolution();
	virtual ~CompactDifferentialEvolution();
	
	double solve(double *xs, int n, int maxEvaluations, int populationSize, double crossover, 
				double mutation, double randomSeed, double (*bounds)[2], double (*F)(double *x));
};

#endif /* COMPACTDIFFERENTIALEVOLUTION_H_ */
