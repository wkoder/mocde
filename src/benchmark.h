/*
 * benchmark.h
 *
 *  Created on: Aug 21, 2011
 *      Author: Moises Osorio
 */

#ifndef BENCHMARK_H_
#define BENCHMARK_H_

#include <string>

using namespace std;

namespace benchmark {
	void setup(char *instanceName, int nreal, int *nobj);
	void evaluate(double *x, double *fx);
	int getEvaluations();
	double **getBounds();
	void destroy();
	double *getVariableDelta();
}

#endif /* BENCHMARK_H_ */
