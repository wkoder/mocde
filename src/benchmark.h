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
	void setup(char *instanceName, int nreal, int *nobj, double **bounds);
	void evaluate(double *x, double *fx);
	int getEvaluations();
	double **getBounds();
}

#endif /* BENCHMARK_H_ */
