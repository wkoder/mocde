/*
 * individual.h
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "config.h"

#if defined MOCDE_IMPL || defined MY_MOEAD_IMPL

#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

class Individual {
public:
	Individual(int nreal, int nobj);
	virtual ~Individual();
	
	double *x;
	double *fx;
	int nreal;
	int nobj;
	
	void copy(Individual *ind);
	Individual *clone();
	
private:

};

#endif /* INDIVIDUAL_H_ */
#endif
