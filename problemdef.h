/*
 * problemdef.h
 *
 *  Created on: October 16, 2011
 *      Author: Moises Osorio
 */

#ifndef PROBLEMDEF_H_
#define PROBLEMDEF_H_

extern int nreal;
extern int nobj;

void wfg1(double *xreal, double *obj);
void wfg2(double *xreal, double *obj);
void wfg6(double *xreal, double *obj);
void dtlz1(double *xreal, double *obj);
void dtlz2(double *xreal, double *obj);
void r_dtlz2(double *xreal, double *obj);
void dtlz3(double *xreal, double *obj);
void dtlz5im(double *xreal, double *obj);
void dtlz7(double *xreal, double *obj);
void zdt1(double *xreal, double *obj);
void zdt2(double *xreal, double *obj);
void zdt3(double *xreal, double *obj);

#endif /* PROBLEMDEF_H_ */
