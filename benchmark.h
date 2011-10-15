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

void benchmarkSetupInstance(int inst, int n, double (*bounds)[2]);
double benchmarkEvaluation(double *x);
int benchmarkGetEvaluations();

#endif /* BENCHMARK_H_ */
