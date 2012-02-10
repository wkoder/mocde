/*
 * util.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#ifndef UTIL_H_
#define UTIL_H_

#include <string>

using namespace std;

#define EPS 1e-9
#define INF 1e20

namespace util {
	string toString(double **x, int n, int m);
	string toString(double *xs, int n);
	double **createMatrix(int r, int c);
	void destroyMatrix(double***m, int r);
}

#endif /* UTIL_H_ */
