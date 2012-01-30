/*
 * util.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <stdio.h>

#include "util.h"

using namespace util;

string util::toString(double **x, int n, int m) {
	string s;
	for (int i = 0; i < n; i++)
		s += toString(x[i], m) + "\n";

	return s;
}

string util::toString(double *xs, int n) {
	string s;
	char buffer[1024];
	for (int i = 0; i < n; i++) {
		if (s.size() > 0)
			s += " ";
		snprintf(buffer, sizeof(buffer), "%.10lf", xs[i]); // Representation of x_i
		s = s.append(buffer);
	}
	
	return s;
}

double **util::createMatrix(int r, int c) {
	double **m = new double*[r];
	for (int i = 0; i < r; i++)
		m[i] = new double[c];
	return m;
}

void util::destroyMatrix(double ***m, int r) {
	for (int i = 0; i < r; i++)
		delete [] (*m)[i];
	delete [] *m;
	*m = NULL;
}
