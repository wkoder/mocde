/*
 * MultiObjectiveCompactDifferentialEvolution.h
 *
 *  Created on: Aug 22, 2011
 *      Author: Moises Osorio
 */

#ifndef MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_
#define MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_

class Individual {
public:
	Individual(int nreal, int nobj);
	virtual ~Individual();
	
	double *x;
	double *fx;
	int nreal;
	int nobj;
	
	void copy(Individual *ind);
	
private:

};

class Subproblem {
public:
	Subproblem(int nreal, int nobj, int populationSize);
	virtual ~Subproblem();
	
	Individual *indiv; // Best solution
	Individual *saved; // Last solution
	double *namda; // Weight vector
	int *table; // Neighbourhood table
private:
	
};

class MultiObjectiveCompactDifferentialEvolution {
public:
	MultiObjectiveCompactDifferentialEvolution();
	virtual ~MultiObjectiveCompactDifferentialEvolution();
	
	int solve(double **xs, double **fxs, int n, int m, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double (*bounds)[2], void (*function)(double *x, double *fx));

private:
	double solve(double *xs, int n, int maxEvaluations, int populationSize, double CR,
					double F, double randomSeed, double (*bounds)[2], double (*function)(double *x));
	
	void initSubproblems(double **L);
	void initNeighborhood();
	void updateIdeal(double *fx);
	int select(int subproblemId, int type);
	void crossover(Individual *a, Individual *b, Individual *c, Individual *child, double CR, double rate);
	void mutation(Individual *ind, double rate);
	void updateSubproblem(Individual *ind, int subproblemId, int type);
	double chebyshevScalarizing(double *fx, double *namda);
	void computeUtility();
	int tourSelection(int *order, int depth);
	
	Subproblem **subproblems;
	int nreal;
	int nobj;
	int populationSize;
	int updateLimit;
	void (*function)(double *x, double *fx);
	double *ideal;
	double *utility;
	int nicheSize;
};

#endif /* MULTIOBJECTIVECOMPACTDIFFERENTIALEVOLUTION_H_ */
