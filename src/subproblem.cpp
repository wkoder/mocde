/*
 * subproblem.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "subproblem.h"

Subproblem::Subproblem(int nreal, int nobj, int nicheSize) {
	namda = new double[nobj];
	table = new int[nicheSize];
	indiv = new Individual(nreal, nobj);
	saved = new Individual(nreal, nobj);
}

Subproblem::~Subproblem() {
	delete [] namda;
	delete [] table;
	delete indiv;
	delete saved;
}
