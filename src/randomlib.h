/*
 * randomlib.h
 *
 *  Created on: Jun 3, 2011
 *      Author: Moises Osorio
 */

#ifndef RANDOMLIB_H_
#define RANDOMLIB_H_

int flip(float prob);
void warmup_random(float random_seed);
int rnd(int low, int high);
int rndint(int n);
float rndreal(float lo, float hi);

#endif /* RANDOMLIB_H_ */
