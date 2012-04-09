// (1+1)-PAES skeleton program code 
/* 
 Copyright (C) 2000, Joshua Knowles and David Corne

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version. 

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details. 

 The GNU General Public License is available at:
 http://www.gnu.org/copyleft/gpl.html
 or by writing to: 
 The Free Software Foundation, Inc., 
 675 Mass Ave, Cambridge, MA 02139, USA.  

 */
//
// PAES is described in several refereed publications. All of them are available 
// for download from Joshua Knowles's webpage: http://dbkgroup.org/knowles/multi/
//
// Please contact the authors if you have any comments, suggestions or questions
// about this file or the Pareto Archived Evolution Strategy (PAES).
// We are at:
// j.knowles@manchester.ac.uk  and  dwcorne@macs.hw.ac.uk 
//
/*********************************************************************************************************
 * To compile : g++ paes.cc -o paes -lm                                                                  *
 * To run :                                                                                              * 
 *  ./paes [problem] [depth] [objectives] [genes] [alleles] [archive] [iterations] [minmax] [pm] [seedpaes]  *
 *                                                                                                       *
 * where all parameters MUST be specified correctly (none are optional) following the instructions below:*
 *                                                                                                       *
 * [problem] - at present there are 3 problems included with this code. They are max, T1, and F5. See    *
 *             below for further details on each of these problems.                                      *
 * [depth] - this is the number of recursive subdivisions of the objective space carried out in order to *
 *           divide the objective space into a grid for the purposes of diversity maintenance. Values of *
 *           between 3 and 6 are useful, depending on number of objectives.                              * 
 * [objectives] - the number of objectives to the problem.                                               *
 * [genes] - the number of genes in the chromosome. (Only integer genes can be accepted).                *
 * [alleles] - the number of values that each gene can take.                                             *
 * [archive] - the maximum number of solutions to be held in the nondominated solutions archive.         *
 * [iterations] - the program terminates after this number of iterations.                                *
 * [minmax] - set this to 0 (zero) for minimization problems, 1 for maximization problems.               *
 * [pm] - the per gene mutation probability. Mutations always change the current allele. Any other       *
 *        viable allele is chosen, with uniform probability.                                             *
 * [seedpaes] - any integer value, to be used as the program's random seedpaes.                                   *
 *                                                                                                       *
 * Example command lines:                                                                                *
 * ./paes max 4 3 20 4 50 10000 1 0.05 3434324   - a valid command line for the max problem              *
 * ./paes T1 5 2 900 2 100 20000 0 0.002 6764421 - a valid command line for the T1 problem               *
 * ./paes F5 4 2 16 16 25 50000 0 0.06 764542    - a valid command line for the F5 problem               *
 *                                                                                                       *
 * The objective function for 3 multiobjective test functions is included in the code.                   *
 *                                                                                                       *
 * The first test function, "max", is a k-objective problem taking a chromosome of k+1 alleles per gene. *
 * It returns, for each objective, the fraction of the chromosome consisting of allele values            *
 * corresponding to the particular objective. i.e. objective 1 counts the number of 1's in the chromosome*
 * and returns this normalised wrt to the length of the chromosome. Similarly for objectives 2, 3, ...   *
 * Zeros in the chromosome do not contribute to any of the objectives. When using the max problem the    *
 * number of alleles must be set to the number of objectives + 1. The problem is a maximization problem  *
 * hence the parameter minmax must be set to 1.                                                          *
 *                                                                                                       *
 * The second objective function, "T1", has been taken from Eckart Zitzler's PhD thesis,                 *
 * Diss. ETH No. 13398 currently available from http://www.tik.ee.ethz.ch/~zitzler/                      *
 * The function is a simple 2-objective minimization                                                     * 
 * problem. The Pareto optimal front is given by 1-sqrt(x) with 1<=x<=0. The code has been written so    *
 * that a chromosome of 900 binary genes is expected, encoding for 30 real numbers between zero and 1.   *
 *                                                                                                       *
 * The third objective function, "F5", was described in "Joshua Knowles and David Corne, The Pareto      *
 * Archived Evolution Strategy: A New Baseline Algortihm for Multiobjective Optimisation, in Proceedings *
 * of the 1999 Congress on Evolutionary Computation (CEC99), pages 98-105, IEEE Service Centre".         *
 * It is a 2-objective problem which accepts a chromosome of k-genes each having k alleles. Briefly the  *
 * function counts the number of adjacent genes having consecutive alleles, reading the chromosome       *
 * forwards (objective 1) and backwards (objective 2). There are k distinct optima in objective space    *
 * but the disribution of these optima in parameter space is highly nonuniform i.e. there are many ways  *
 * of obtaining some of the optima (they are easy to find) and very few ways of finding others.          *
 * The problem is cast as a minimization problem. For a chromosome of 16 genes with 16 alleles the       *
 * problem is already very difficult.                                                                    *
 *                                                                                                       *
 * In the code presented here, maintenance of diversity is carried out in the objective space. i.e.      *
 * solutions which are genetically different but which give rise to the same objective values are not    *
 * allowed to co-exist in the archive. The grid used for maintaining diversity also takes into account   *
 * objective values only. Parameter-space, and genotypic measures of diversity may also be used with     *
 * PAES but they have not been included in this skeleton code.                                           *              
 *                                                                                                       *
 ********************************************************************************************************/

#include "paes.h"

#ifdef PAES_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "../benchmark.h"

#define MAX_GENES 1000  // change as necessary
#define MAX_OBJ 10     // change as necessary 
#define MAX_POP 10     // change as necessary
#define MAX_ARC 200    // change as necessary
#define MAX_LOC 32768  // number of locations in grid (set for a three-objective problem using depth 5)
#define LARGE 2000000000 // should be about the maximum size of an integer for your compiler 
FILE *fp;

typedef struct solution {
	int chrom[MAX_GENES];
	double obj[MAX_OBJ];
	int grid_loc;
} sol;

sol *c; // current solution
sol *m; // mutant solution
sol *arc; // archive of solutions

int compare_min(double *first, double *second, int n);
int compare_max(double *first, double *second, int n);
int equal(double *first, double *second, int n);
void init();
void getX(sol *s, double *x);
void print_genome(sol *s);
void print_eval(sol *s);
void evaluate(sol *s);
void max(sol *s);
void T1(sol *s);
void F5(sol *s);
void UF1(sol *s);
void add_to_archive(sol *s);
void mutate(sol *s);
int compare_to_archive(sol *s);
void update_grid(sol *s);
void archive_soln(sol *s);
int find_loc(double *eval);

double pm; // mutation rate as a decimal probability per bit (zero probability of generating current allele, uniform probability of generating all other allele values. i.e. for binary chromosomes this is the "flip" mutation rate.) 
int seedpaes, depth, genes, alleles, archive, objectives, iterations, minmax; // command line parameters 
int arclength = 0; // current length of the archive

double gl_offset[MAX_OBJ]; // the range, offset etc of the grid 
double gl_range[MAX_OBJ];
double gl_largest[MAX_OBJ];
int grid_pop[MAX_LOC]; // the array holding the population residing in each grid location

void (*problem)(double *x, double *fx);
double **bounds;
int parameters;

int paes(double **xb, double **fxb, double **_bounds, void (*_problem)(double *x, double *fx), int _depth, int _parameters, int _objectives, 
		int _genes, int _alleles, int _archive, int _iterations, int _minmax, double _pm, int _seed) {
	int i;
	int result;

	bounds = _bounds;
	problem = _problem;
	depth =_depth;
	parameters = _parameters;
	objectives = _objectives;
	genes = _genes;
	alleles = _alleles;
	archive = _archive;
	iterations = _iterations;
	minmax = _minmax;
	pm = _pm;
	seedpaes = _seed;

	srand(seedpaes); // seedpaes system random number generator

	// begin (1+1)-PAES
	init();
	//  print_genome(c);    // Uncomment to check genome is correct 
	evaluate(c);
	//  print_eval(c);      // Uncomment to check objective values generated
	//  getchar();  
	add_to_archive(c);
	// begin main loop
	for (i = 0; i < iterations-1; i++) {
		//if (i%100==0)            // just print out the number of iterations done every 100
		//printf("%d\n", i);
		*m = *c; // copy the current solution 
		mutate(m); // and mutate using the per-bit mutation rate specified by the command param pm

		//   print_genome(m);
		evaluate(m);

		//print_eval(m);
		// print_genome(m);
		//  printf("Comparing ");
		//print_eval(c); 
		// printf("and ");
		//print_eval(m);

		//MINIMIZE MAXIMIZE
		if (minmax == 0)
			result = compare_min(c->obj, m->obj, objectives);
		else
			result = compare_max(c->obj, m->obj, objectives);

		// printf("RESULT = %d\n", result);
		// printf("arclength = %d\n", arclength);
		if (result != 1) { // if mutant is not dominated by current (else discard it)
			if (result == -1) { // if mutant dominates current
				//printf("m dominates c\n");
				update_grid(m); //calculate grid location of mutant solution and renormalize archive if necessary
				archive_soln(m); //update the archive by removing all dominated individuals
				*c = *m; // replace c with m
			} else if (result == 0) { // if mutant and current are nondominated wrt each other
				result = compare_to_archive(m);
				if (result != -1) { // if mutant is not dominated by archive (else discard it)
					update_grid(m);
					archive_soln(m);
					// if mutant dominates the archive or is in less crowded grid loc than c
					if ((grid_pop[m->grid_loc] <= grid_pop[c->grid_loc]) || (result == 1)) {  
						*c = *m; // then replace c with m
					}
				}
			}
		}
		/*  printf("\nArchive = \n");
		 for (j = 0; j < arclength; j++)
		 print_eval(&arc[j]);
		 getchar();*/

	}
	//printf("\nThe Archive is now... \n");
	for (i = 0; i < arclength; i++) {
//		print_eval(&arc[i]);
		getX(&arc[i], xb[i]);
		for (int j = 0; j < objectives; j++)
			fxb[i][j] = arc[i].obj[j];
	}

	/*  printf("\nThe genetic material in the archive is now... \n");
	 for (i = 0; i < arclength; i++)
	 print_genome(&arc[i]);
	 */
	
	return arclength;
}

int compare_to_archive(sol *s) // compares a solution to every member of the archive. Returns -1 if dominated by
		{ // any member, 1 if dominates any member, and 0 otherwise
	int i = 0;
	int result = 0;

	while ((i < arclength) && (result != 1) && (result != -1)) {
		//MINIMIZE MAXIMIZE
		if (minmax == 0)
			result = compare_min(m->obj, (&arc[i])->obj, objectives);
		else
			result = compare_max(m->obj, (&arc[i])->obj, objectives);
		i++;
	}

	return (result);
}

void init() {
	int i;
	int val;

	// allocate memory for solutions
	c = (sol *) malloc(MAX_POP * sizeof(sol));
	m = (sol *) malloc(MAX_POP * sizeof(sol));
	arc = (sol *) malloc(MAX_ARC * sizeof(sol));

	if ((!c) || (!m) || (!arc)) {
		printf("Out of memory. Aborting.\n");
		exit(-1);
	}

	// initialise c with a uniform distribution of values from 0 to alleles-1
	c->grid_loc = 0;
	for (i = 0; i < genes; i++) {
		val = (int) (alleles * (rand() / (RAND_MAX + 1.0)));
		c->chrom[i] = val;
	}
}

void add_to_archive(sol *s) {
	arc[arclength] = *s;
	arclength++;
}

void mutate(sol *s) // apply mutation to chromosome s using per-gene mutation probability pm.  
		{
	int i;
	int var;

	for (i = 0; i < genes; i++) {
		if ((rand() / (RAND_MAX + 1.0)) < pm) // mutate gene ?
				{
			var = 1 + (int) ((alleles - 1) * (rand() / (RAND_MAX + 1.0))); // generate var between 1 and alleles-2
			s->chrom[i] = (s->chrom[i] + var) % alleles; // add var to allele value of current gene and mod.  
		}
	}
}

void print_genome(sol *s) {
	int i;
	for (i = 0; i < genes; i++)
		printf("%d ", s->chrom[i]);
	printf("\n");
}

void print_eval(sol *s) {
	int i;
	for (i = 0; i < objectives; i++)
		printf("%g ", s->obj[i]);
	printf("\n");
}

void getX(sol *s, double *x) {
	if (genes % parameters != 0) {
		printf("Genes is not a multiple of parameters! Exiting...\n");
		exit(0);
	}
	
	int genesByParam = genes / parameters;
	for (int i = 0; i < parameters; i++) {
		double val = 0;
		for (int j = 0; j < genesByParam; j++)
			val = 2*val + s->chrom[i*genesByParam + j];
		
		val /= pow(2.0, genesByParam) - 1;
		x[i] = bounds[i][0] + val*(bounds[i][1] - bounds[i][0]);
	}
}

void evaluate(sol *s) {
	double x[parameters];
	
	getX(s, x);
	problem(x, s->obj);
}

void max(sol *s) {
	int i;
	double sc[MAX_OBJ + 1];

	if (alleles != objectives + 1) {
		printf("You've given invalid command line parameters for function max. Check paes.cc for details. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < objectives; i++)
		sc[i] = 0;

	for (i = 0; i < genes; i++)
		sc[s->chrom[i]]++;

	for (i = 1; i <= objectives; i++) {
		sc[i] /= (double) genes;
		s->obj[i - 1] = sc[i];
	}
}

void T1(sol *s) {
	// test function T1 described in Eckart Zitzler's PhD thesis. See notes at top for further details 
	int i, j;
	int mul;
	double var[30];
	double sum = 0;
	double f1, f2, g, h;

	if ((genes != 900) || (alleles != 2) || (objectives != 2)) {
		printf(
				"T1 requires 900 binary genes. You have %d genes and they are %d-ary.\nIt is a 2 objective problem. You have %d objectives. Exiting.\n",
				genes, alleles, objectives);
		exit(-1);
	}

	//set paramaters to zero
	for (i = 0; i < 30; i++)
		var[i] = 0;

	//convert from binary into integer for the 30 params
	for (i = 0; i < 30; i++) {
		mul = 1;
		for (j = 29; j >= 0; j--) {
			var[i] += mul * s->chrom[i * 30 + j];
			mul *= 2;
		}
	}

	// normalize the params between 0 and 1
	for (i = 0; i < 30; i++)
		var[i] /= 1 << 30;

	f1 = var[0];

	for (i = 1; i < 30; i++)
		sum += var[i];
	g = 1 + (9 * sum) / 29.0;

	h = 1 - sqrt(f1 / g);

	f2 = g * h;

	s->obj[0] = f1;
	s->obj[1] = f2;

}

void F5(sol *s) // two-objective minimization problem
		{
	int i;

	if ((objectives != 2) || (genes != alleles)) {
		printf("You've given invalid command parameters for function F5. Check paes.cc for details. Exiting.\n");
		exit(-1);
	}

	s->obj[0] = genes - 1; // worst score for each objective is genes-1
	s->obj[1] = genes - 1;

	for (i = 0; i < genes - 1; i++) {
		if (s->chrom[i + 1] == s->chrom[i] + 1) // reduce score of objective 1 if there are adjacent genes having consecutive values
			s->obj[0]--;

		if (s->chrom[i + 1] == s->chrom[i] - 1) // as above but reading in reverse for objective 2
			s->obj[1]--;
	}
}

void UF1(sol *s) {
	int j, count1, count2;
	double sum1, sum2, yj;

	sum1 = sum2 = 0.0;
	count1 = count2 = 0;
	for (j = 2; j <= genes; j++) {
		yj = 1.0 * s->chrom[j - 1] / alleles - sin(6.0 * M_PI * s->chrom[0] / alleles + j * M_PI / genes);
		yj = yj * yj;
		if (j % 2 == 0) {
			sum2 += yj;
			count2++;
		} else {
			sum1 += yj;
			count1++;
		}
	}
	s->obj[0] = 1.0 * s->chrom[0] / alleles + 2.0 * sum1 / (double) count1;
	s->obj[1] = 1.0 - sqrt(1.0 * s->chrom[0] / alleles) + 2.0 * sum2 / (double) count2;
}

int compare_min(double *first, double *second, int n) {
	// compares two n-dimensional vectors of objective values for minimization problems
	// returns 1 if first dominates second,
	// -1 if second dominates first, and 0 otherwise

	int obj = 0;
	int deflt = 0;
	int current;

	do {

		if (*first < *second)
			current = 1;
		else if (*second < *first)
			current = -1;
		else
			current = 0;
		if ((current) && (current == -deflt)) {

			return (0);
		}
		if (current != 0) {
			deflt = current;
		}
		obj++;
		*first++;
		*second++;
	} while (obj < n);

	return (deflt);
}

int compare_max(double *first, double *second, int n) {
	// as for compare_min() but for maximization problems 

	int obj = 0;
	int deflt = 0;
	int current;

	do {
		if (*first > *second)
			current = 1;
		else if (*second > *first)
			current = -1;
		else
			current = 0;
		if ((current) && (current == -deflt)) {
			return (0);
		}
		if (current != 0) {
			deflt = current;
		}
		obj++;
		*first++;
		*second++;
	} while (obj < n);
	return (deflt);
}

int equal(double *first, double *second, int n) {
	// checks to n-dimensional vectors of objectives to see if they are identical
	// returns 1 if they are, 0 otherwise

	int obj = 0;

	do {
		if (*first != *second)
			return (0);
		*first++;
		*second++;
		obj++;
		// printf("%d\n",obj);
	} while (obj < n);
	return (1);
}

void archive_soln(sol *s) {
	// given a solution s, add it to the archive if
	// a) the archive is empty 
	// b) the archive is not full and s is not dominated or equal to anything currently in the archive
	// c) s dominates anything in the archive                        
	// d) the archive is full but s is nondominated and is in a no more crowded square than at least one solution
	// in addition, maintain the archive such that all solutions are nondominated. 

	int i;
	int repl = 0;
	int yes = 0;
	int most;
	int result;
	int join = 0;
	int old_arclength;
	int set = 0;

	int tag[MAX_ARC];
	sol *tmp;
	if (!(tmp = (sol *) malloc(MAX_ARC * sizeof(sol)))) {
		printf("Out of memory\n");
		exit(-1);
	}

	for (i = 0; i < archive; i++) {
		tag[i] = 0;
	}

	if (arclength == 0) {
		add_to_archive(s);
		return;
	}

	i = 0;
	result = 0;
	while ((i < arclength) && (result != -1)) {
		result = equal(s->obj, (&arc[i])->obj, objectives);
		if (result == 1)
			break;
		//MINIMIZE MAXIMIZE
		if (minmax == 0)
			result = compare_min(s->obj, (&arc[i])->obj, objectives);
		else
			result = compare_max(s->obj, (&arc[i])->obj, objectives);
		//  printf("%d\n", result);

		if ((result == 1) && (join == 0)) {
			arc[i] = *s;
			join = 1;
		} else if (result == 1) {
			tag[i] = 1;
			set = 1;
		}
		i++;
	}

	old_arclength = arclength;
	if (set == 1) {
		for (i = 0; i < arclength; i++) {
			tmp[i] = arc[i];
		}
		arclength = 0;

		for (i = 0; i < old_arclength; i++) {
			if (tag[i] != 1) {
				arc[arclength] = tmp[i];
				arclength++;
			}
		}
	}

	if ((join == 0) && (result == 0)) { // ie solution is non-dominated by the list
		if (arclength == archive) {
			most = grid_pop[s->grid_loc];
			for (i = 0; i < arclength; i++) {
				if (grid_pop[(&arc[i])->grid_loc] > most) {
					most = grid_pop[(&arc[i])->grid_loc];
					repl = i;
					yes = 1;
					//   printf("i = %d\n", i);
				}
			}
			if (yes) {
				arc[repl] = *s;
			}
		} else {
			add_to_archive(s);
		}
	}
	free(tmp);
}

int find_loc(double *eval) {
	// find the grid location of a solution given a vector of its objective values

	int loc = 0;
	int d;
	int n = 1;

	int i;

	int inc[MAX_OBJ];
	double width[MAX_OBJ];

	// printf("obj = %d, depth = %d\n", objectives, depth);

	// if the solution is out of range on any objective, return 1 more than the maximum possible grid location number
	for (i = 0; i < objectives; i++) {
		if ((eval[i] < gl_offset[i]) || (eval[i] > gl_offset[i] + gl_range[i]))
			return ((int) pow(2.0, (objectives * depth)));
	}

	for (i = 0; i < objectives; i++) {
		inc[i] = n;
		n *= 2;
		width[i] = gl_range[i];
	}

	for (d = 1; d <= depth; d++) {
		for (i = 0; i < objectives; i++) {
			if (eval[i] < width[i] / 2 + gl_offset[i])
				loc += inc[i];
			else
				gl_offset[i] += width[i] / 2;
		}
		for (i = 0; i < objectives; i++) {
			inc[i] *= (objectives * 2);
			width[i] /= 2;
		}
	}
	return (loc);
}

void update_grid(sol *s) {
	// recalculate ranges for grid in the light of a new solution s
	static int change = 0;
	int a, b;
	int square;
	double offset[MAX_OBJ];
	double largest[MAX_OBJ];
	double sse;
	double product;

	for (a = 0; a < objectives; a++) {
		offset[a] = LARGE;
		largest[a] = -LARGE;
	}

	for (b = 0; b < objectives; b++) {
		for (a = 0; a < arclength; a++) {
			if ((&arc[a])->obj[b] < offset[b])
				offset[b] = (&arc[a])->obj[b];
			if ((&arc[a])->obj[b] > largest[b])
				largest[b] = (&arc[a])->obj[b];
		}
	}
//    printf("oldCURENT:largest = %f, offset = %f\n", largest[0], offset[0]); 
//    printf("oldCURENT:largest = %f, offset = %f\n", largest[1], offset[1]); 

	for (b = 0; b < objectives; b++) {
		if (s->obj[b] < offset[b])
			offset[b] = s->obj[b];
		if (s->obj[b] > largest[b])
			largest[b] = s->obj[b];
	}

	sse = 0;
	product = 1;

	for (a = 0; a < objectives; a++) {

		sse += ((gl_offset[a] - offset[a]) * (gl_offset[a] - offset[a]));
		sse += ((gl_largest[a] - largest[a]) * (gl_largest[a] - largest[a]));
		product *= gl_range[a];
	}

	// printf("sse = %f\n", sse);
	if (sse > (0.1 * product * product)) //if the summed squared error (difference) between old and new
										 //minima and maxima in each of the objectives
			{ //is bigger than 10 percent of the square of the size of the space
		change++; // then renormalise the space and recalculte grid locations

		for (a = 0; a < objectives; a++) {
			gl_largest[a] = largest[a] + 0.2 * largest[a];
			gl_offset[a] = offset[a] + 0.2 * offset[a];
			gl_range[a] = gl_largest[a] - gl_offset[a];
		}

		for (a = 0; a < pow(2.0, (objectives * depth)); a++) {
			grid_pop[a] = 0;
		}

		for (a = 0; a < arclength; a++) {
			square = find_loc((&arc[a])->obj);
//			printf("SQUARE: %d\n", square);
			(&arc[a])->grid_loc = square;
			grid_pop[square]++;

		}
	}
	square = find_loc(s->obj);
//	printf("SQUARE: %d\n", square);
	s->grid_loc = square;
	grid_pop[(int) pow(2.0, (objectives * depth))] = -5;
	grid_pop[square]++;

}

#endif
