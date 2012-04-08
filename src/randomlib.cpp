/*
 * randomlib.c
 */

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/types.h>

#include "randomlib.h"

static double oldrand[55];                               /* Array of 55 random numbers */
static int jrand;                                             /* current random number */
static double rndx2;                                /* used with random normal deviate */
static int rndcalcflag;                             /* used with random normal deviate */

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

/* Get seed number for random and start it up */
void randomize(double seed)
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random (seed);
    return;
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

/* Rutina para generacion de numeros aleatorios gaussianos con
   media cero y desviacion estandar sigma.

   Codigo original de Joachim Sprave, University of Dortmund, 
   Alemania, Dept. of Computer Science and Systems Analysis

*/

#define NRAND_SAMPLES 5
#define Uniform randreal

static int Seed = 0;
static int isinit = 0;

double Gauss(double sigma)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double u, x, y, u0, u1, u2;


/* 	SIGMA	--> standard deviation */

/* L1: */
    u = Uniform();
    u0 = Uniform();
    if (u >= .919544406) {
	goto L2;
    }
    x = (u0 + u * .825339283) * 2.40375766 - 2.11402808;
    goto L10;
  L2:
    if (u < .965487131) {
	goto L4;
    }
  L3:
    u1 = Uniform();
    y = sqrt(4.46911474 - log(u1) * 2.);
    u2 = Uniform();
    if (y * u2 > 2.11402808) {
	goto L3;
    }
    goto L9;
  L4:
    if (u < .949990709) {
	goto L6;
    }
  L5:
    u1 = Uniform();
    y = u1 * .273629336 + 1.84039875;
    u2 = Uniform();
    if (exp(y * -.5 * y) * .39894228 - .443299126 + y * .209694057 < u2 *
	.0427025816) {
	goto L5;
    }
    goto L9;
  L6:
    if (u < .925852334) {
	goto L8;
    }
  L7:
    u1 = Uniform();
    y = u1 * 1.55066917 + .289729574;
    u2 = Uniform();
    if (exp(y * -.5 * y) * .39894228 - .443299126 + y * .209694057 < u2 *
	.0159745227) {
	goto L7;
    }
    goto L9;
  L8:
    u1 = Uniform();
    y = u1 * .289729574;
    u2 = Uniform();
    if (exp(y * -.5 * y) * .39894228 - .382544556 < u2 * .0163977244) {
	goto L8;
    }
  L9:
    x = y;
    if (u0 >= .5) {
	x = -y;
    }
  L10:
    ret_val = sigma * x;
    return ret_val;
}

double normal(double m, double sigma)
{
    return m + Gauss(sigma);
}

void initrandom(int seed)
{
    Seed = seed;
}

double randreal(void)
{
    double result;

    if (!isinit) {
	srand(Seed);
	isinit = 1;
    }
    result = ((double) rand());
    result /= RAND_MAX;

    return (result);
}

int randint(int lo, int hi)
{
    return (lo + randreal() * (hi - lo + 1));
}
