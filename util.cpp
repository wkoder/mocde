/*
 * util.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <stdio.h>

#include "util.h"

string toString(double *xs, int n) {
	string s;
	char buffer[1024];
	for (int i = 0; i < n; i++) {
		if (s.size() > 0)
			s += ", ";
		snprintf(buffer, sizeof(buffer), "%.10lf", xs[i]); // Representation of x_i
		s = s.append(buffer);
	}
	
	return s;
}
