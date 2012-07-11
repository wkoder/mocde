/*
 * MultiObjectiveCompactDifferentialEvolution.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#include "mocde.h"

#ifdef MOCDE_IMPL

#include <stdio.h>
#include <math.h>
#include <boost/math/special_functions/erf.hpp>

#include "randomlib.h"
#include "benchmark.h"
#include "stats.h"
#include "util.h"
#include "moead/global.h"

using namespace util;

static double SQRT_2 = sqrt(2.0);

MultiObjectiveCompactDifferentialEvolution::MultiObjectiveCompactDifferentialEvolution() {
	
}

MultiObjectiveCompactDifferentialEvolution::~MultiObjectiveCompactDifferentialEvolution() {

}

double MultiObjectiveCompactDifferentialEvolution::sampleValue(double mu, double sigma) {
	double sqrt2Sigma = SQRT_2 * sigma;
	double erfMuNeg = boost::math::erf((mu - 1) / sqrt2Sigma);
	double erfMuPlus = boost::math::erf((mu + 1) / sqrt2Sigma);
	
	double u = randreal();
	double C = - boost::math::erf((mu + 1) / sqrt2Sigma) / (erfMuNeg - erfMuPlus);
	return mu - sqrt2Sigma * boost::math::erf_inv((u - C) * (erfMuNeg - erfMuPlus));
}

void MultiObjectiveCompactDifferentialEvolution::denormalizeSolution(double *x, double *normx, double **bounds, int nreal) {
	for (int i = 0; i < nreal; i++) {
		double len = bounds[i][1] - bounds[i][0];
		x[i] = bounds[i][0] + len/2 + normx[i]*len/2;
	}
}

void MultiObjectiveCompactDifferentialEvolution::sampleSolution(double *normx, double *u, double *d, int nreal) {
	for (int i = 0; i < nreal; i++) {
		normx[i] = sampleValue(u[i], d[i]);
		if (normx[i] < -1 || normx[i] > 1) {
			cerr << "ERROR WITH SAMPLING! " << normx[i] << " = u: " << u[i] << " d: " << d[i] << endl;
		}
	}
}

double MultiObjectiveCompactDifferentialEvolution::getCrowdDistance(Individual *ind, vector<Individual *> &archive) {
	double indDist = INF;
	for (vector<Individual *>::iterator it = archive.begin(); it < archive.end(); it++)
		if (*it != ind)
			indDist = min(indDist, norm2((*it)->fx, ind->fx, nobj));
	
	return indDist;
}

bool MultiObjectiveCompactDifferentialEvolution::addToArchive(vector<Individual *> &archive, Individual *ind, Individual *parent) {
	if (true)
		return addToMOEADArchive(archive, ind, parent);
	
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
				double F, int maxSurvival, double randomSeed, double **bounds, void (*function)(double *x, double *fx)) {
	this->nreal = nreal;
	this->nobj = nobj;
	this->populationSize = populationSize;
	this->function = function;
	
	double u[nreal];
	double d[nreal];
	double xr[nreal];
	double xs[nreal];
	double xt[nreal];
	vector<Individual *> archive;
	
	initMOEADArchive(archive);
	
	// PV initialization
	for (int i = 0; i < nreal; i++) {
		u[i] = 0;
		d[i] = 10;
	}

	Individual *off = new Individual(nreal, nobj);
	double offNormX[nreal];
	
	Individual *elite = new Individual(nreal, nobj);
	double eliteNormX[nreal];
	sampleSolution(eliteNormX, u, d, nreal);
	denormalizeSolution(elite->x, eliteNormX, bounds, nreal);
	function(elite->x, elite->fx);
	addToArchive(archive, elite, elite);
	int survivedIterations = 0; // Elitism level, persistent when > maxEvaluations
	while (benchmark::getEvaluations() < maxEvaluations) {
//		cout << "U: " << util::toString(u, nreal) << endl;
//		cout << "D: " << util::toString(d, nreal) << endl;
//		cout << "E: " << util::toString(eliteNormX, nreal) << endl << endl;
		stats::report(archive);
		
//		if (benchmark::getEvaluations() % 1000 == 0) {
//			cerr << "Ideal: ";
//			for (int i = 0; i < nobj; i++)
//				cerr << ideal[i] << " ";
//			cerr << endl;
//			cerr << "Dev: ";
//			for (int i = 0; i < nreal; i++)
//				cerr << d[i] << " ";
//			cerr << endl;
//		}
		
		// Mutation
		// TODO Get sampling stats
#ifdef RESAMPLING
		do {
#endif
		sampleSolution(xr, u, d, nreal);
		sampleSolution(xs, u, d, nreal);
		sampleSolution(xt, u, d, nreal);
		
		bool wasUnfeasible = false;
		for (int i = 0; i < nreal; i++) {
#ifdef RAND_BEST_1
			offNormX[i] = xt[i] + F*(xr[i] - xs[i]) + F*(eliteNormX[i] - xt[i]);
#else
			offNormX[i] = xt[i] + F*(xr[i] - xs[i]);
#endif
			if (offNormX[i] < -1) {
				offNormX[i] = -1;
				wasUnfeasible = true;
			} else if (offNormX[i] > 1) {
				offNormX[i] = 1;
				wasUnfeasible = true;
			}
		}
#ifdef RESAMPLING
		while (wasUnfeasible);
#endif
		
		// Crossover
		for (int i = 0; i < nreal; i++)
			if (!flip(CR))
				offNormX[i] = eliteNormX[i];
		
		// Elite selection
		denormalizeSolution(off->x, offNormX, bounds, nreal);
		function(off->x, off->fx);
		
		double *winner = eliteNormX;
		double *loser = offNormX;
		ParetoDominance cmp = comparePareto(off->fx, elite->fx, nobj);
		bool replaceElite = false;
		// TODO: Test Pareto dominance removal
		if (cmp == DOMINATES || survivedIterations >= maxSurvival) {
			addToArchive(archive, off, elite);
			winner = offNormX;
			loser = eliteNormX;
			replaceElite = true;
			survivedIterations = 0;
		} else if (cmp == NON_DOMINATED && addToArchive(archive, off, elite)) {
			winner = offNormX;
			loser = eliteNormX;
			replaceElite = true;
			survivedIterations = 0;
		} else
			survivedIterations++;
			
		// PV update TODO: Check if dominance is needed!
		for (int i = 0; i < nreal; i++) {
			double u2 = u[i] + (winner[i] - loser[i]) / populationSize;
			u2 = max(u2, -1.0);
			u2 = min(u2, 1.0);
			double d2 = sqrt(fabs(d[i]*d[i] + u[i]*u[i] - u2*u2 + 
					(winner[i]*winner[i] - loser[i]*loser[i]) / populationSize));
			
			u[i] = u2;
			d[i] = d2;
		}
		
		if (replaceElite) {
			elite->copy(off);
			for (int i = 0; i < nreal; i++)
				eliteNormX[i] = offNormX[i];
		}
	}
	stats::report(archive);
	
	delete elite;
	delete off;
	
	util::removeDominated(archive);
	util::getPopulationX(xb, archive);
	util::getPopulationFx(fxb, archive);
	return archive.size();
}

void MultiObjectiveCompactDifferentialEvolution::initMOEADArchive(std::vector<Individual *> &archive) {
	ideal = new double[nobj];
	L = util::createMatrix(populationSize, nobj);
	
	for (int i = 0; i < nobj; i++)
		ideal[i] = INF;
	
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
}

bool MultiObjectiveCompactDifferentialEvolution::addToMOEADArchive(std::vector<Individual *> &archive, Individual *ind, Individual *parent) {
	updateIdeal(ind->fx);
	if (archive.size() == 0) {
		for (int i = 0; i < populationSize; i++)
			archive.push_back(ind->clone());
		
		return true;
	}
	
	bool added = false;
	for (int i = 0; i < populationSize; i++) {
		Individual *sub = archive[i];
		double f1 = chebyshevScalarizing(sub->fx, L[i]);
		double f2 = chebyshevScalarizing(ind->fx, L[i]);
//		cerr << f1 << " " << f2 << endl;
		if (f2 < f1) {
			added = true;
			sub->copy(ind);
		}
	}
	
	return added;
}

double MultiObjectiveCompactDifferentialEvolution::chebyshevScalarizing(double *fx, double *namda) {
	double max = -INF;

	for (int n = 0; n < nobj; n++) {
		double diff = fabs(fx[n] - ideal[n]);
		double feval = diff * (namda[n] < EPS ? 0.0001 : namda[n]);
		if (max < feval)
			max = feval;
	}

	return max;
}

void MultiObjectiveCompactDifferentialEvolution::updateIdeal(double *fx) {
	for (int i = 0; i < nobj; i++)
		if (fx[i] < ideal[i])
			ideal[i] = fx[i];
}

#endif
