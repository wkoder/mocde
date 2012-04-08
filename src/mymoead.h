/*
 * Moead.h
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "config.h"

#ifdef MY_MOEAD_IMPL

#ifndef MY_MOEAD_H_
#define MY_MOEAD_H_

#include "subproblem.h"
#include "individual.h"

class MyMOEAD {
public:
	MyMOEAD();
	virtual ~MyMOEAD();
	int solve(double **xb, double **fxb, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
					double F, double (*bounds)[2]);

private:
	Subproblem **subproblems;
	int updateLimit;
	double *ideal;
	double *utility;
	int nicheSize;
	int nobj;
	int nreal;
	int populationSize;
	
	double chebyshevScalarizing(double *fx, double *namda);
	void initSubproblems(double **L);
	void initNeighborhood();
	void updateIdeal(double *fx);
	int select(int subproblemId, int type);
	void crossover(Individual *a, Individual *b, Individual *c, Individual *child, double CR, double rate);
	void mutation(Individual *ind, double rate);
	void updateSubproblem(Individual *ind, int subproblemId, int type);
	void computeUtility();
	int tourSelection(int *order, int depth);
};

#endif /* MY_MOEAD_H_ */
#endif
