/*
 * PAES_main.h
 *
 *  Created on: Apr 8, 2012
 *      Author: Moises Osorio
 */

#include "../config.h"

#ifdef JMETAL

#ifndef PAES_MAIN_H_
#define PAES_MAIN_H_

class PAES_main {
public:
	
	PAES_main();
	virtual ~PAES_main();
	
	int solve(double **xs, double **fxs, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
					double F, double randomSeed, double **bounds);
};


#endif /* PAES_MAIN_H_ */
#endif
