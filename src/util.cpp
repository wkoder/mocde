/*
 * util.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "util.h"

namespace util {

string toString(double **x, int n, int m) {
	string s;
	for (int i = 0; i < n; i++)
		s += toString(x[i], m) + "\n";

	return s;
}

string toString(int x) {
	string s;
	char buffer[20];
	snprintf(buffer, sizeof(buffer), "%d", x);
	s = s.append(buffer);
	return s;
}

string toString(double *xs, int n) {
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

double **createMatrix(int r, int c) {
	double **m = new double*[r];
	for (int i = 0; i < r; i++)
		m[i] = new double[c];
	return m;
}

void destroyMatrix(double ***m, int r) {
	for (int i = 0; i < r; i++)
		delete [] (*m)[i];
	delete [] *m;
	*m = NULL;
}

void printx(const char *d, double *x, int n) {
	printf("%s = ", d);
	for (int i = 0; i < n; i++)
		printf("%.3f ", x[i]);
	printf("\n");
}

double norm2(double *a, double *b, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	return sum;
}

bool comparePair(pair<double, int> a, pair<double, int> b) {
	return a.first < b.first;
}

ParetoDominance comparePareto(double *a, double *b, int nobj) {
	ParetoDominance cmp = EQUAL;
	for (int i = 0; i < nobj; i++)
		if (fabs(a[i] - b[i]) < EPS)
			continue;
		else if (a[i] < b[i]) // Dominates on f_i
			if (cmp == DOMINATED)
				return NON_DOMINATED;
			else
				cmp = DOMINATES;
		else //if (a[i] > b[i]) // Dominated on f_i
			if (cmp == DOMINATES)
				return NON_DOMINATED;
			else
				cmp = DOMINATED;
	
	return cmp;
}

void removeDominated(std::vector<Individual *> &population) {
	for (unsigned int i = 0; i < population.size(); i++)
		for (unsigned int j = i+1; j < population.size(); j++) {
			ParetoDominance cmp = comparePareto(population[i]->fx, population[j]->fx, population[i]->nobj);
			if (cmp == DOMINATED) {
				population.erase(population.begin() + i);
				i--;
				break;
			}
		}
}

void writeMatrixFile(double **matrix, int rows, int cols, string filename) {
	ofstream psfile(filename.c_str());
	if (psfile.is_open()) {
		psfile << toString(matrix, rows, cols);
		psfile.close();
	} else
		cout << "Cannot open file " << filename.c_str() << endl;
}

void getPopulationX(double **x, vector<Individual *> population) {
	for (unsigned int i = 0; i < population.size(); i++) {
		Individual *ind = population[i];
		copyVector(x[i], ind->x, ind->nreal);
	}
}

void getPopulationFx(double **fx, vector<Individual *> population) {
	for (unsigned int i = 0; i < population.size(); i++) {
		Individual *ind = population[i];
		copyVector(fx[i], ind->fx, ind->nobj);
	}
}

void copyVector(double *copy, double *x, int len) {
	for (int i = 0; i < len; i++)
		copy[i] = x[i];
}

}
