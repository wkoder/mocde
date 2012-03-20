/*
 * MultiObjectiveCompactDifferentialEvolution.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <math.h>
#include <stdio.h>

#include "mocde.h"
#include "randomlib.h"
#include "benchmark.h"
#include "util.h"

#include "moead/algorithm.h"
#include "paes/paes.h"
#include "nsga2/global.h"
#include "mymoead.h"

using namespace util;

#define MOEAD_IMPL
//#define MY_MOEAD_IMPL
//#define MOCDE_IMPL
//#define PAES_IMPL
//#define NSGA2_IMPL

MultiObjectiveCompactDifferentialEvolution::MultiObjectiveCompactDifferentialEvolution() {
	
}

MultiObjectiveCompactDifferentialEvolution::~MultiObjectiveCompactDifferentialEvolution() {

}

void ensureBounds(double *x, double (*bounds)[2], int n) {
	for (int i = 0; i < n; i++) {
		x[i] = max(x[i], bounds[i][0]);
		x[i] = min(x[i], bounds[i][1]);
	}
}

void generateX(double *x, double *u, double *d, double (*bounds)[2], int n) {
	for (int i = 0; i < n; i++) {
		x[i] = normal(u[i], d[i]);
		while (x[i] < bounds[i][0] || x[i] > bounds[i][1]) // Until it's tight
			x[i] = normal(u[i], d[i]);
	}
//	ensureBounds(x, bounds, n);
}

double getCrowdDistance(Individual *ind, vector<Individual *> &archive) {
	double indDist = INF;
	for (vector<Individual *>::iterator it = archive.begin(); it < archive.end(); it++)
		if (*it != ind)
			indDist = min(indDist, norm2((*it)->fx, ind->fx, nobj));
	
	return indDist;
}

/**
 * Implementation of the Fisher-Yates shuffling algorithm 
 * (http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle)
 */
void shuffle(int *x, int n) {
	for (int i = n-1; i >= 0; i--) {
		int j = rndint(i+1);
		int tmp = x[i];
		x[i] = x[j];
		x[j] = tmp;
	}
}

bool MultiObjectiveCompactDifferentialEvolution::addToArchive(vector<Individual *> &archive, Individual *ind, Individual *parent) {
	bool dominates = false;
	vector<Individual *>::iterator it = archive.begin();
	while (it < archive.end()) {
		ParetoDominance cmp = comparePareto(ind->fx, (*it)->fx, nobj);
		if (cmp == DOMINATED)
			return false;
		
		if (cmp == DOMINATES) {
			delete *it;
			it = archive.erase(it);
			dominates = true;
		} else
			it++;
	}
	
	if (dominates) {
		archive.push_back(ind->clone());
		return true;
	}
	
	if (archive.size() < (unsigned int) populationSize) {
		archive.push_back(ind->clone());
	} else {
		double indDist = getCrowdDistance(ind, archive);
		double crowdedDist = INF;
		Individual *crowdedInd = NULL;
		for (it = archive.begin(); it < archive.end(); it++) {
			double dist = getCrowdDistance(*it, archive);
			if (dist < crowdedDist) {
				crowdedDist = dist;
				crowdedInd = *it;
			}
		}
		
		if (crowdedDist < indDist) { // Improves diversity
			archive.erase(find(archive.begin(), archive.end(), crowdedInd));
			delete crowdedInd;
			archive.push_back(ind->clone());
		}
	}
	
	return getCrowdDistance(ind, archive) > getCrowdDistance(parent, archive); // Accept if less crowded than parent
}


int MultiObjectiveCompactDifferentialEvolution::solve(double **xb, double **fxb, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx)) {
	this->nreal = nreal;
	this->nobj = nobj;
	this->populationSize = populationSize;
	this->function = function;
	
	randomize(randomSeed);
	initrandom(randomSeed * (1 << 30));
	
//#ifdef MOCDE_IMPL
//	int effortForIdeal = populationSize/10;
//	int budget = maxEvaluations / (populationSize + nobj*(effortForIdeal-1));
//	int newPopulationSize = 100;
//	double u[nreal];
//	double d[nreal];
//	for (int i = nobj-1; i >= 0; i--) {
//		for (int j = 0; j < nreal; j++) {
//			double r = bounds[j][1] - bounds[j][0];
//			u[j] = bounds[j][0] + r/2;
//			d[j] = r * 0.341;
//		}
//		
//		currentFullNamdaIdx = i;
//		ideal[i] = solve(xb[i], nreal, budget*effortForIdeal, newPopulationSize, CR, F, u, d, bounds, chebyshevSimpleCost);
//		cout << ideal[i] << " ";
//	}
//	cout << endl;
//	
//	double dx = 1.0 / (populationSize - 1);
//	for (int i = 0; i < populationSize; i++) {
//		L[i] = new double[2];
//		L[i][0] = i * dx;
//		L[i][1] = 1 - L[i][0];
//	}
//	
//	for (int sub = 2; sub < populationSize; sub++) {
//		for (int j = 0; j < nreal; j++) {
//			double r = bounds[j][1] - bounds[j][0];
//			u[j] = bounds[j][0] + r/2;
//			d[j] = r * 0.341;
//		}
//		
//		currentNamda = L[sub-1];
//		solve(xb[sub], nreal, budget, newPopulationSize, CR, F, u, d, bounds, chebyshevCost);
//	}
//	
//	for (int i = 0; i < populationSize; i++) // Get the real function value
//		function(xb[i], fxb[i]);
//#endif
#ifdef MOCDE_IMPL
	double u[nreal];
	double d[nreal];
	double xr[nreal];
	double xs[nreal];
	double xt[nreal];
	vector<Individual *> archive;
	
	// PV initialization
	for (int i = 0; i < nreal; i++) {
		double r = bounds[i][1] - bounds[i][0];
		u[i] = bounds[i][0] + r/2;
		d[i] = r*4;
	}
	
	Individual *ind = new Individual(nreal, nobj);
	generateX(ind->x, u, d, bounds, nreal);
	function(ind->x, ind->fx);
	archive.push_back(ind->clone());
	Individual *off = new Individual(nreal, nobj);
	while (benchmark::getEvaluations() < maxEvaluations) {
		// Mutation
		generateX(xr, u, d, bounds, nreal);
		generateX(xs, u, d, bounds, nreal);
		generateX(xt, u, d, bounds, nreal);
		
		for (int i = 0; i < nreal; i++) {
			off->x[i] = xt[i] + F*(xr[i] - xs[i]);
			// Ensure bounds
			if (off->x[i] < bounds[i][0])
				off->x[i] = bounds[i][0];
			else if (off->x[i] > bounds[i][1])
				off->x[i] = bounds[i][1];
//			if (off->x[i] < bounds[i][0])
//				off->x[i] = rndreal(bounds[i][0], xt[i]);
//			else if (off->x[i] > bounds[i][1])
//				off->x[i] = rndreal(xt[i], bounds[i][1]);
		}
		
		// Crossover
//		Individual *elite = archive[rndint(archive.size())];
		for (int i = 0; i < nreal; i++)
			if (!flip(CR))
				off->x[i] = ind->x[i];
		
		// Elite selection
		function(off->x, off->fx);
		Individual *winner = ind;
		Individual *loser = off;
		ParetoDominance cmp = comparePareto(off->fx, ind->fx, nobj);
		if (cmp == DOMINATES) {
			addToArchive(archive, off, ind);
			winner = off;
			loser = ind;
		} else if (cmp == NON_DOMINATED) {
			if (addToArchive(archive, off, ind)) {
				winner = off;
				loser = ind;
			}
		}
			
		// PV update
		for (int i = 0; i < nreal; i++) {
			double u2 = u[i] + (winner->x[i] - loser->x[i]) / populationSize;
			// TODO: Improve bound ensuring!
			u2 = max(u2, bounds[i][0]);
			u2 = min(u2, bounds[i][1]);

			double d2 = sqrt(fabs(d[i]*d[i] + u[i]*u[i] - u2*u2 + 
					(winner->x[i]*winner->x[i] - loser->x[i]*loser->x[i]) / populationSize));
			
			u[i] = u2;
//			if (d2 > d[i]+EPS)
//				cout << "YES! " << d2-d[i] << endl;
			d[i] = d2;
		}
		
		if (winner == off)
			ind->copy(off);
	}
	
//	printx("Mean", u, nreal);
//	printx("Dev", d, nreal);
	
	delete ind;
	delete off;
	
	for (unsigned int i = 0; i < archive.size(); i++) {
		for (int j = 0; j < nreal; j++)
			xb[i][j] = archive[i]->x[j];
		for (int j = 0; j < nobj; j++)
			fxb[i][j] = archive[i]->fx[j];
	}
	
	return archive.size();
#endif
#ifdef PAES_IMPL
	return paes(xb, fxb, bounds, function, 4, nreal, nobj, nreal*30, 2, populationSize, maxEvaluations, 0, 0.01, (int)(randomSeed*(1<<30)));
#endif
#ifdef NSGA2_IMPL
	return nsga2(xb, fxb, populationSize, maxEvaluations/populationSize, 0.9, 0.01, nreal, nobj, bounds, randomSeed);
#endif
#ifdef MOEAD_IMPL
//	seed = 177;
	seed = (int)(randomSeed * (1 << 30));
	rnd_uni_init = 90.0;
	lowBound = 0;
	uppBound = 1;
	double **L = createMatrix(populationSize, nobj);
	int nicheSize = nobj == 2 ? 100 : 150;
	int updateLimit = nicheSize / 10;
		
	char filename[1024];
	sprintf(filename, "resources/W%dD.dat", nobj);
	ifstream file(filename);
	if (!file.is_open()) {
		cout << "File " << filename << " not found.\n";
		exit(-1);
	}
	for (int i = 0; i < populationSize; i++)
		for (int j = 0; j < nobj; j++)
			file >> L[i][j];
	file.close();
	
	CMOEAD MOEAD;
	MOEAD.load_parameter(populationSize, maxEvaluations/populationSize, nicheSize, updateLimit, F, maxEvaluations);
	return MOEAD.exec_emo(xb, fxb, L);
#endif
#ifdef MY_MOEAD_IMPL
	MyMOEAD moead;
	return moead.solve(xb, fxb, nobj, maxEvaluations, populationSize, CR, F, randomSeed, bounds);
#endif
}
