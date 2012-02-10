/*
 * randomlib.h
 *
 *  Created on: Jun 3, 2011
 *      Author: Moises Osorio
 */

#ifndef RANDOMLIB_H_
#define RANDOMLIB_H_

float randomperc();
int flip(float prob);
void randomize(double seed);
void warmup_random(float random_seed);
int rnd(int low, int high);
int rndint(int n);
float rndreal(float lo, float hi);

double Gauss(double sigma);
double normal(double m, double sigma);
void initrandom(int seed);
double randreal(void);
int randint(int lo, int hi);

#endif /* RANDOMLIB_H_ */
