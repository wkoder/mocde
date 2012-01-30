/*
 * paes.h
 *
 *  Created on: Dec 5, 2011
 *      Author: Moises Osorio [WCoder]
 */

#ifndef PAES_H_
#define PAES_H_

int paes(double **xb, double **fxb, double (*_bounds)[2], void (*_problem)(double *x, double *fx), int _depth, int _parameters, int _objectives, 
		int _genes, int _alleles, int _archive, int _iterations, int _minmax, double _pm, int _seed);

#endif /* PAES_H_ */
