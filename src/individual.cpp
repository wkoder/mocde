/*
 * individual.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "individual.h"

Individual::Individual(int nreal, int nobj) {
	this->nreal = nreal;
	this->nobj = nobj;
	this->x = new double[nreal];
	this->fx = new double[nobj];
}

Individual::~Individual() {
	delete [] x;
	delete [] fx;
}

Individual *Individual::clone() {
	Individual *clone = new Individual(nreal, nobj);
	clone->copy(this);
	
	return clone;
}

void Individual::copy(Individual *ind) {
	for (int i = 0; i < nreal; i++)
		x[i] = ind->x[i];
	for (int i = 0; i < nobj; i++)
		fx[i] = ind->fx[i];
}
