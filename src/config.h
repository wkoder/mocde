/*
 * config.h
 * 
 * Configuration file to know which implementation is going to be compiled.
 *
 *  Created on: Apr 7, 2012
 *      Author: Moises Osorio
 */

#ifndef CONFIG_H_
#define CONFIG_H_

//#define MOEAD_IMPL
//#define MY_MOEAD_IMPL
#define MOCDE_IMPL
//#define MOCDE_SCALAR
//#define PAES_IMPL
//#define NSGA2_IMPL

//#define TOOLKIT_FUN // Functions (mainly WFG) from the Toolkit directory

//#define RESAMPLING // If boundary constraint is being handled by re-sampling the individual, otherwise projection is used
#define RAND_BEST_1
//#define RAND_BEST_2
//#define RAND_1
//#define OUTPUT_INTERVAL 5000
#define OUTPUT_INTERVAL 1e9

#ifdef PAES_IMPL
#define JMETAL // Implementations from JMetal: http://jmetal.sourceforge.net
#endif

#endif /* CONFIG_H_ */
