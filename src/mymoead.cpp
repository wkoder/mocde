/*
 * Moead.cpp
 *
 *  Created on: Feb 10, 2012
 *      Author: Moises Osorio [WCoder]
 */

#include "mymoead.h"

#ifdef MY_MOEAD_IMPL

#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "util.h"
#include "subproblem.h"
#include "randomlib.h"
#include "benchmark.h"

using namespace std;
using namespace util;

MyMOEAD::MyMOEAD() {

}

MyMOEAD::~MyMOEAD() {
	
}

//MyMOEAD *currentInstance;
//
//double *currentNamda;
//int currentFullNamdaIdx;
//
//double chebyshevSimpleCost(double *x) {
//	double fx[currentInstance->nobj];
//	currentInstance->function(x, fx);
//	
//	double sum = 0;
//	for (int i = 0; i < currentInstance->nobj; i++)
//		if (i == currentInstance->nobj)
//			sum += fx[i];
//		else
//			sum += fx[i] * 0.0001;
//	return sum;
//}
//
//double chebyshevCost(double *x) {
//	double fx[currentInstance->nobj];
//	currentInstance->function(x, fx);
//	return currentInstance->chebyshevScalarizing(fx, currentNamda);
//}

double MyMOEAD::chebyshevScalarizing(double *fx, double *namda) {
	double max = -INF;

	for (int n = 0; n < nobj; n++) {
		double diff = fabs(fx[n] - ideal[n]);
		double feval = diff * (namda[n] < EPS ? 0.0001 : namda[n]);
		if (max < feval)
			max = feval;
	}

	return max;
}

void MyMOEAD::updateIdeal(double *fx) {
	for (int i = 0; i < nobj; i++)
		if (ideal[i] > fx[i])
			ideal[i] = fx[i];
}

void MyMOEAD::initSubproblems(double **L) {
	for (int i = 0; i < nobj; i++)
		ideal[i] = INF;
	for (int i = 0; i < populationSize; i++)
		utility[i] = 1.0;
	
	subproblems = new Subproblem*[populationSize];
	for (int i = 0; i < populationSize; i++) {
		Subproblem *sub = new Subproblem(nreal, nobj, nicheSize);
		for (int j = 0; j < nreal; j++) // Random initialization
			sub->indiv->x[j] = rndreal(bounds[j][0], bounds[j][1]);
		for (int j = 0; j < nreal; j++)
			if (sub->indiv->x[j] < -5 || sub->indiv->x[j] > 5)
				cerr << "FUCKKKKKKK! " << j << ": " << sub->indiv->x[j] << " not in [" << bounds[j][0] << ", " << bounds[j][1] << "]" << endl;
		benchmark::evaluate(sub->indiv->x, sub->indiv->fx); // Evaluate individual
		updateIdeal(sub->indiv->fx); // Update reference point
		sub->saved->copy(sub->indiv);
		for (int j = 0; j < nobj; j++)
			sub->namda[j] = L[i][j];
		subproblems[i] = sub;
	}
}

void MyMOEAD::initNeighborhood() {
	pair<double, int> dist[populationSize];
	for (int i = 0; i < populationSize; i++) {
		for (int j = 0; j < populationSize; j++) {
			dist[j].first = norm2(subproblems[i]->namda, subproblems[j]->namda, nicheSize);
			dist[j].second = j;
		}
		
		sort(dist, dist+populationSize, comparePair);
		for (int j = 0; j < nicheSize; j++) // TODO: Avoid adding itself
			subproblems[i]->table[j] = dist[j].second;
	}
}

int MyMOEAD::select(int subproblemId, int type) {
	return type == 1 ? subproblems[subproblemId]->table[rndint(nicheSize)] : rndint(populationSize);
}

void MyMOEAD::updateSubproblem(Individual *ind, int subproblemId, int type) {
	int time = 0;
	int size = type == 1 ? nicheSize : populationSize; // From neighborhood or from the whole population
//	int perm[size];
//	for (int i = 0; i < size; i++)
//		perm[i] = i;
//	shuffle(perm, size); // Shuffle the permutation - WTF? Not needed, I think it was replaced by tourSelection()

	for (int i = 0; i < size; i++) {
		// Pick a subproblem to update
//		int k = type == 1 ? subproblems[subproblemId]->table[perm[i]] : perm[i];
		int k = type == 1 ? subproblems[subproblemId]->table[i] : i;

		// Calculate the values of objective function regarding the current subproblem
		double f1 = chebyshevScalarizing(subproblems[k]->indiv->fx, subproblems[k]->namda);
		double f2 = chebyshevScalarizing(ind->fx, subproblems[k]->namda);
//		cout << f2 << " < " << f1 << " = " << (f2<fMyMOEAD::solve1) << endl;
		if (f2 < f1) {
			subproblems[k]->indiv->copy(ind);
			time++;
		}
		
		// The maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= updateLimit)
			break;
	}
}

void MyMOEAD::computeUtility() {
	for (int n = 0; n < populationSize; n++) {
		double f1 = chebyshevScalarizing(subproblems[n]->indiv->fx, subproblems[n]->namda);
		double f2 = chebyshevScalarizing(subproblems[n]->saved->fx, subproblems[n]->namda);
		if (f1 > f2) {
//			cout << "ERROR! f1: " << f1 << ", f2: " << f2 << endl;
//			exits(-1);
		}
		
		double delta = (f2 - f1) / f2;
//		cout << utility[n] << " -> ";
		if (delta > 0.001) {
			utility[n] = 1.0;
		} else {
			//uti = 0.95*(1.0 + delta/0.001) * utility[n];
			//utility[n] = uti < 1.0 ? uti : 1.0;
			utility[n] = (0.95 + 0.05 * delta / 0.001) * utility[n];
		}
		
		subproblems[n]->saved->copy(subproblems[n]->indiv);
//		cout << utility[n] << " (" << delta << "), ";
	}
//	cout << endl;
}

int MyMOEAD::tourSelection(int *order, int depth) {
	// Selection based on utility
	int orderSize;
	for (orderSize = 0; orderSize < nobj; orderSize++)
		order[orderSize] = orderSize; // Select first m weights

	vector<int> candidate;
	for (int n = nobj; n < populationSize; n++)
		candidate.push_back(n); // Set of unselected weights

	for (; orderSize < floor(populationSize / 5.0); orderSize++) {
		int best_idd = rndint(candidate.size());
		int best_sub = candidate[best_idd];
		for (int i = 1; i < depth; i++) {
			int i2 = rndint(candidate.size());
			int s2 = candidate[i2];
			if (utility[s2] > utility[best_sub]) {
				best_idd = i2;
				best_sub = s2;
			}
		}
		
		order[orderSize] = best_sub;
		candidate.erase(candidate.begin() + best_idd);
	}
	
	return orderSize;
}

void MyMOEAD::crossover(Individual *a, Individual *b, Individual *c, Individual *child, double CR, double rate) {
	for (int i = 0; i < nreal; i++) {
		if (flip(CR)) // Crossover
			child->x[i] = a->x[i] + rate*(c->x[i] - b->x[i]);
		else
			child->x[i] = a->x[i];
		
		if (child->x[i] < bounds[i][0])
			child->x[i] = rndreal(bounds[i][0], a->x[i]);
		if (child->x[i] > bounds[i][1])
			child->x[i] = rndreal(a->x[i], bounds[i][1]);
	}
}

void MyMOEAD::mutation(Individual *ind, double rate) {
	double eta_m = 20;
	for (int j = 0; j < nreal; j++)
		if (flip(rate)) {
			double y = ind->x[j];
			double yl = bounds[j][0];
			double yu = bounds[j][1];
			double delta1 = (y - yl) / (yu - yl);
			double delta2 = (yu - y) / (yu - yl);
			double rnd = rndreal(0, 1);
			double mut_pow = 1.0 / (eta_m + 1.0);
			double deltaq;
			if (rnd <= 0.5) {
				double xy = 1.0 - delta1;
				double val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			} else {
				double xy = 1.0 - delta2;
				double val = 2.0 * (1.0 - rnd) + 2.0 * (rnd - 0.5) * (pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			
			y = y + deltaq * (yu - yl);
			if (y < yl)
				y = yl;
			if (y > yu)
				y = yu;
			ind->x[j] = y;
		}
}

int MyMOEAD::solve(double **xb, double **fxb, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
				double F, double **bounds) {
	this->nobj = nobj;
	this->nreal = nreal;
	this->populationSize = populationSize;
	this->ideal = new double[nobj];
	this->nicheSize = nobj == 2 ? 100 : 150;
	this->updateLimit = this->nicheSize / 10;
	this->utility = new double[populationSize];
	this->bounds = bounds;
	double **L = createMatrix(populationSize, nobj);
	
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
	
	initSubproblems(L);
	initNeighborhood();
	
	int order[populationSize];
	int generation = 0;
	while (benchmark::getEvaluations() < maxEvaluations) {
		int orderSize = tourSelection(order, 10);
		for (int i = 0; i < orderSize; i++) {
			int subId = order[i];
			int type = flip(F) ? 1 : 2;
			int a = select(subId, type);
			while (a == subId)
				a = select(subId, type);
			int b = select(subId, type);
			while (a == b || subId == b)
				b = select(subId, type);
			
			// Create child
			Individual *child = new Individual(nreal, nobj);
			crossover(subproblems[subId]->indiv, subproblems[a]->indiv, subproblems[b]->indiv, child, 1.0, 0.5);
			mutation(child, 1.0 / nreal); // Mutate child
			
			benchmark::evaluate(child->x, child->fx); // Evaluate child
			updateIdeal(child->fx); // Update reference point
			updateSubproblem(child, subId, type); // Update subproblem
			
			delete child;
		}
		
//			for (int i = 0; i < nobj; i++)
//				cout << ideal[i] << " ";
//			cout << endl;
		
		generation++;
		if (generation % 50 == 0) // Every 50 generations
			computeUtility();
	}
	
	for (int i = 0; i < populationSize; i++) {
		for (int j = 0; j < nreal; j++)
			xb[i][j] = subproblems[i]->indiv->x[j];
		for (int j = 0; j < nobj; j++)
			fxb[i][j] = subproblems[i]->indiv->fx[j];
		
		delete subproblems[i];
	}
	
	delete [] subproblems;
	delete [] utility;
	delete [] ideal;
	destroyMatrix(&L, populationSize);
	
	return populationSize;
}

#endif
