/*
 * util.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <vector>
#include <string>
#include "individual.h"

using namespace std;

#define EPS 1e-9
#define INF 1e20

enum ParetoDominance {
	DOMINATED, 
	DOMINATES, 
	EQUAL, 
	NON_DOMINATED 
};

namespace util {
	string toString(double **x, int n, int m);
	string toString(double *xs, int n);
	double **createMatrix(int r, int c);
	void destroyMatrix(double***m, int r);
	void printx(const char *d, double *x, int n);
	double norm2(double *a, double *b, int n);
	bool comparePair(pair<double, int> a, pair<double, int> b);
	ParetoDominance comparePareto(double *a, double *b, int nobj);
	void removeDominated(vector<Individual *> &population);
}

#endif /* UTIL_H_ */
