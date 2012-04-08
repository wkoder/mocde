/*
 * subproblem.h
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "config.h"

#ifdef MY_MOEAD_IMPL

#ifndef SUBPROBLEM_H_
#define SUBPROBLEM_H_

#include "individual.h"

class Subproblem {
public:
	Subproblem(int nreal, int nobj, int populationSize);
	virtual ~Subproblem();
	
	Individual *indiv; // Best solution
	Individual *saved; // Last solution
	double *namda; // Weight vector
	int *table; // Neighbourhood table
	
private:
	
};

#endif /* SUBPROBLEM_H_ */
#endif
