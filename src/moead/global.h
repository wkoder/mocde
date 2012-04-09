#include "../config.h"

#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>

#include "random.h"
#include "../benchmark.h"

using namespace std;

//------------- Parameters in test instance ------------------

extern int nreal;
extern int nobj;                    //  the number of variables and objectives

extern double lowBound;
extern double uppBound;   //  lower and upper bounds of variables

//extern char    strTestInstance[256];




//------------- Parameters in random number ------------------
extern int     seed;
extern long    rnd_uni_init;        


//------------- Parameters in MOEA/D -------------------------

extern vector <double> idealpoint;
extern double          scale[100];  


extern int etax;
extern int etam;   // distribution indexes of crossover and mutation

extern double real;
extern double realm;
extern double realb;   // crossover, mutation, selection probabilities


#endif
