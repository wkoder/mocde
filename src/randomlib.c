/*
 * randomlib.c
 */

#include <math.h>

#include "randomlib.h"

//float  Rseed;			       		         /* Random numbers seed */
double oldrand[55];                               /* Array of 55 random numbers */
int jrand;                                             /* current random number */
double rndx2;                                /* used with random normal deviate */
int rndcalcflag;                             /* used with random normal deviate */

/* Create next batch of 55 random numbers */
void advance_random() {
	int j1;
	double new_random;
	
	for (j1 = 0; j1 < 24; j1++) {
		new_random = oldrand[j1] - oldrand[j1 + 31];
		if (new_random < 0.0)
			new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
	
	for (j1 = 24; j1 < 55; j1++) {
		new_random = oldrand[j1] - oldrand[j1 - 24];
		if (new_random < 0.0)
			new_random = new_random + 1.0;
		oldrand[j1] = new_random;
	}
	
}

/* initialization routine for randomnormaldeviate */
void initrandomnormaldeviate() {
	rndcalcflag = 1;
}

/* normal noise with specified mean & std dev: mu & sigma */
double noise(double mu, double sigma) {
	double randomnormaldeviate();
	return ((randomnormaldeviate() * sigma) + mu);
}

/* Get random off and running */
void warmup_random(float random_seed) {
	int j1, ii;
	double new_random, prev_random;
	
	oldrand[54] = random_seed;
	new_random = 0.000000001;
	prev_random = random_seed;
	
	for (j1 = 1; j1 <= 54; j1++) {
		ii = (21 * j1) % 54;
		oldrand[ii] = new_random;
		new_random = prev_random - new_random;
		if (new_random < 0.0)
			new_random = new_random + 1.0;
		prev_random = oldrand[ii];
	}

	advance_random();
	advance_random();
	advance_random();
	
	jrand = 0;
}

///* Initialize random numbers batch */
//void randomize() {
//	int j1;
//	
//	for (j1 = 0; j1 <= 54; j1++)
//		oldrand[j1] = 0.0;
//	
//	jrand = 0;
//	warmup_random(Rseed);
//}

/* Fetch a single random number between 0.0 and 1.0 - Subtractive Method */
/* See Knuth, D. (1969), v. 2 for details */
/* name changed from random() to avoid library conflicts on some machines*/
float randomperc() {
	jrand++;
	
	if (jrand >= 55) {
		jrand = 1;
		advance_random();
	}
	
	return ((float) oldrand[jrand]);
}

/* random normal deviate after ACM algorithm 267 / Box-Muller Method */
double randomnormaldeviate() {
	float randomperc();
	double t, rndx1;
	
	if (rndcalcflag) {
		rndx1 = sqrt(-2.0 * log((double) randomperc()));
		t = 6.2831853072 * (double) randomperc();
		rndx2 = sin(t);
		rndcalcflag = 0;
		return (rndx1 * cos(t));
	} else {
		rndcalcflag = 1;
		return (rndx2);
	}
}

/* Pick a random integer between low and high */
int rnd(int low, int high) {
	int i;
	float randomperc();
	
	if (low >= high)
		i = low;
	else {
		i = (randomperc() * (high - low + 1)) + low;
		if (i > high)
			i = high;
	}
	
	return (i);
}

/* real random number between specified limits */
float rndreal(float lo, float hi) {
	return ((randomperc() * (hi - lo)) + lo);
}

/* Picks a random integer between 0 and n-1 (inclusive) */
int rndint(int n) {
	return rnd(0, n-1);
}

/* Flip a biased coin - true if heads */
int flip(float prob) {
	if (randomperc() <= prob)
		return (1);
	else
		return (0);
	
}
