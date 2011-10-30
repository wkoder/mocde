/*
 * MultiObjectiveCompactDifferentialEvolution.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include <math.h>
#include <stdio.h>
#include <iostream>

#include "mocde.h"
#include "rand.h"
#include "randomlib.h"
#include "benchmark.h"
#include "util.h"

#include "moead/algorithm.h"

#undef EPS
#define EPS 1e-9

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
	for (int i = 0; i < n; i++)
		x[i] = normal(u[i], d[i]);
	ensureBounds(x, bounds, n);
}

void printx(char *d, double *x, int n) {
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

double MultiObjectiveCompactDifferentialEvolution::chebyshevScalarizing(double *fx, double *namda) {
	double max = -1e30;

	for (int n = 0; n < nobj; n++) {
		double diff = fabs(fx[n] - ideal[n]);
		double feval = diff * (namda[n] < EPS ? 0.0001 : namda[n]);
		if (max < feval)
			max = feval;
	}

	return max;
}

double MultiObjectiveCompactDifferentialEvolution::solve(double *xb, int n, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], double (*function)(double *x)) {
	double u[n];
	double d[n];
	double elite[n];
	double xr[n];
	double xs[n];
	double xt[n];
	double xoff[n];
	warmup_random(randomSeed);
	initrandom(randomSeed * (1 << 30));
	int startingEval = benchmark::getEvaluations();
	
	// PV initialization
	for (int i = 0; i < n; i++) {
		double r = bounds[i][1] - bounds[i][0];
		u[i] = bounds[i][0] + r/2;
		d[i] = r * 0.341;
	}
	
	generateX(elite, u, d, bounds, n);
	double felite = function(elite);
	while (benchmark::getEvaluations()-startingEval < maxEvaluations) {
//		if (benchmark::getEvaluations() % populationSize == 0) {
//			printf("Iteration #%d:\n", benchmark::getEvaluations() / populationSize);
//			printf("	best: %s\n", util::toString(elite, n).c_str());
//			printf("	f(best): %.6f\n", felite);
//		}
		
		// Mutation
		generateX(xr, u, d, bounds, n);
		generateX(xs, u, d, bounds, n);
		generateX(xt, u, d, bounds, n);
		
		for (int i = 0; i < n; i++)
			xoff[i] = xt[i] + F*(xr[i] - xs[i]);
		ensureBounds(xoff, bounds, n);
		
		// Crossover
		for (int i = 0; i < n; i++)
			if (!flip(CR))
				xoff[i] = elite[i];
		
		// Elite selection
		double fxoff = function(xoff);
		double *winner = elite;
		double *loser = xoff;
		if (fxoff < felite) {
			winner = xoff;
			loser = elite;
			felite = fxoff;
		}
		
		// PV update
		for (int i = 0; i < n; i++) {
			double u2 = u[i] + (winner[i] - loser[i]) / populationSize;
			double d2 = sqrt(fabs(d[i]*d[i] + u[i]*u[i] - u2*u2 + 
					(winner[i]*winner[i] - loser[i]*loser[i]) / populationSize));
			
			u2 = max(u2, bounds[i][0]);
			u2 = min(u2, bounds[i][1]);
			u[i] = u2;
			d[i] = d2;
		}
		
		if (winner == xoff)
			for (int i = 0; i < n; i++)
				elite[i] = xoff[i];
	}
	
//	printf("Iteration #%d:\n", benchmark::getEvaluations() / populationSize);
//	printf("	best: %s\n", util::toString(elite, n).c_str());
//	printf("	f(best): %.6f\n", felite);
	
	for (int i = 0; i < n; i++)
		xb[i] = elite[i];
	
	return felite;
}

void MultiObjectiveCompactDifferentialEvolution::updateIdeal(double *fx) {
	for (int i = 0; i < nobj; i++)
		if (ideal[i] > fx[i])
			ideal[i] = fx[i];
}

void MultiObjectiveCompactDifferentialEvolution::initSubproblems(double **L) {
	for (int i = 0; i < nobj; i++)
		ideal[i] = 1e30;
	for (int i = 0; i < populationSize; i++)
		utility[i] = 1.0;
	
	subproblems = new Subproblem*[populationSize];
	for (int i = 0; i < populationSize; i++) {
		Subproblem *sub = new Subproblem(nreal, nobj, nicheSize);
		for (int j = 0; j < nreal; j++) // Random initialization
			sub->indiv->x[j] = rndreal(0, 1);
		function(sub->indiv->x, sub->indiv->fx); // Evaluate individual
		updateIdeal(sub->indiv->fx); // Update reference point
		sub->saved->copy(sub->indiv);
		for (int j = 0; j < nobj; j++)
			sub->namda[j] = L[i][j];
		subproblems[i] = sub;
	}
}

void MultiObjectiveCompactDifferentialEvolution::initNeighborhood() {
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

int MultiObjectiveCompactDifferentialEvolution::select(int subproblemId, int type) {
	return type == 1 ? subproblems[subproblemId]->table[rndint(nicheSize)] : rndint(populationSize);
}

void MultiObjectiveCompactDifferentialEvolution::crossover(Individual *a, Individual *b, Individual *c, Individual *child, double CR, double rate) {
	for (int i = 0; i < nreal; i++) {
		if (flip(CR)) // Crossover
			child->x[i] = a->x[i] + rate*(c->x[i] - b->x[i]);
		else
			child->x[i] = a->x[i];
		
//		if (child->x[i] < xbounds[i][0])
//			child->x[i] = rndreal(xbounds[i][0], a->x[i]);
//		if (child->x[i] > xbounds[i][1])
//			child->x[i] = r
		if (child->x[i] < 0)
			child->x[i] = rndreal(0, a->x[i]);
		if (child->x[i] > 1)
			child->x[i] = rndreal(a->x[i], 1);
	}
}

void MultiObjectiveCompactDifferentialEvolution::mutation(Individual *ind, double rate) {
	double eta_m = 20;
	for (int j = 0; j < nreal; j++)
		if (flip(rate)) {
			double y = ind->x[j];
			double yl = 0;
			double yu = 1;
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

void MultiObjectiveCompactDifferentialEvolution::updateSubproblem(Individual *ind, int subproblemId, int type) {
	int time = 0;
	int size = type == 1 ? nicheSize : populationSize; // From neighborhood or from the whole population
	int perm[size];
	for (int i = 0; i < size; i++)
		perm[i] = i;
	shuffle(perm, size); // Shuffle the permutation

	for (int i = 0; i < size; i++) {
		// Pick a subproblem to update
		int k = type == 1 ? subproblems[subproblemId]->table[perm[i]] : perm[i];

		// Calculate the values of objective function regarding the current subproblem
		double f1 = chebyshevScalarizing(subproblems[k]->indiv->fx, subproblems[k]->namda);
		double f2 = chebyshevScalarizing(ind->fx, subproblems[k]->namda);
//		cout << f2 << " < " << f1 << " = " << (f2<fMultiObjectiveCompactDifferentialEvolution::solve1) << "\n";
		if (f2 < f1) {
			subproblems[k]->indiv->copy(ind);
			time++;
		}
		
		// The maximal number of solutions updated is not allowed to exceed 'limit'
		if (time >= updateLimit)
			break;
	}
}

void MultiObjectiveCompactDifferentialEvolution::computeUtility() {
	for (int n = 0; n < populationSize; n++) {
		double f1 = chebyshevScalarizing(subproblems[n]->indiv->fx, subproblems[n]->namda);
		double f2 = chebyshevScalarizing(subproblems[n]->saved->fx, subproblems[n]->namda);
		if (f1 > f2) {
//			cout << "ERROR! f1: " << f1 << ", f2: " << f2 << "\n";
//			exits(-1);
		}
		double delta = (f2 - f1) / f2;
		//double delta = (f2 - f1);
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
//	cout << "\n";
}

int MultiObjectiveCompactDifferentialEvolution::tourSelection(int *order, int depth) {
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

int MultiObjectiveCompactDifferentialEvolution::solve(double **xb, double **fxb, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx)) {
	this->nreal = nreal;
	this->nobj = nobj;
	this->populationSize = populationSize;
	this->function = function;
	this->ideal = new double[nobj];
	this->nicheSize = nobj == 2 ? 100 : 150;
	this->updateLimit = this->nicheSize / 10;
	this->utility = new double[populationSize];
	
	warmup_random(randomSeed);
	initrandom(randomSeed * (1 << 30));
	
	double **L = util::createMatrix(populationSize, nobj);
	char filename[1024];
	sprintf(filename, "moead/W%dD.dat", nobj);
	
	ifstream file(filename);
	if (!file.is_open()) {
		cout << "File " << filename << " not found.\n";
		exit(-1);
	}
	for (int i = 0; i < populationSize; i++)
		for (int j = 0; j < nobj; j++)
			file >> L[i][j];
	file.close();
//			double dx = 1.0 / (populationSize - 1);
//			for (int i = 0; i < populationSize; i++) {
//					L[i] = new double[2];
//					if (i % 2 == 0)
//						L[i][0] = i/2 * dx;
//					else
//						L[i][0] = L[i-1][1];
//					L[i][1] = 1 - L[i][0];
//			}
	
	
	if (false) {
		seed = 177;
		rnd_uni_init = 90.0;
		lowBound = 0;
		uppBound = 1;
		
		CMOEAD MOEAD;
		MOEAD.load_parameter(populationSize, maxEvaluations/populationSize, nicheSize, updateLimit, F);
		MOEAD.exec_emo(xb, fxb, L);
	} else {
		initSubproblems(L);
		initNeighborhood();
		
		int order[populationSize];
		int generation = 0;
		while (benchmark::getEvaluations() < maxEvaluations) {
//			printx("x", subproblems[0]->indiv->x, nreal);
//			printx("xs", subproblems[0]->saved->x, nreal);
//			cout << "\n";
			
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
				
				function(child->x, child->fx); // Evaluate child
				updateIdeal(child->fx); // Update reference point
				updateSubproblem(child, subId, type); // Update subproblem
				
				delete child;
			}
			
//			for (int i = 0; i < nobj; i++)
//				cout << ideal[i] << " ";
//			cout << "\n";
			
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
	}
	
	util::destroyMatrix(&L, populationSize);
	return populationSize;
}

Individual::Individual(int nreal, int nobj) {
	this->nreal = nreal;
	this->nobj = nobj;
	this->x = new double[nreal];
	this->fx = new double[nobj];
}

Individual::~Individual() {
	delete [] x;
	delete [] fx;
}

void Individual::copy(Individual *ind) {
	for (int i = 0; i < nreal; i++)
		x[i] = ind->x[i];
	for (int i = 0; i < nobj; i++)
		fx[i] = ind->fx[i];
}

Subproblem::Subproblem(int nreal, int nobj, int nicheSize) {
	namda = new double[nobj];
	table = new int[nicheSize];
	indiv = new Individual(nreal, nobj);
	saved = new Individual(nreal, nobj);
}

Subproblem::~Subproblem() {
	delete [] namda;
	delete [] table;
	delete indiv;
	delete saved;
}
