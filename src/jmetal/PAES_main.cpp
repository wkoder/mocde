#include "PAES_main.h"

#ifdef JMETAL

//  PAES_main.cpp
//
//  Author:
//       Cristian Zambrano <cristian_uteq@hotmail.com>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "BenchmarkProblem.h"
#include "Algorithm.h"
#include "Operator.h"
#include "Solution.h"
#include "RealSolutionType.h"
#include "PolynomialMutation.h"
#include <iostream>
#include "paesj.h"
#include "Real.h"
//#include <time.h>

PAES_main::PAES_main() {
	
}

PAES_main::~PAES_main() {
	
}

int PAES_main::solve(double **xs, double **fxs, int nreal, int nobj, int maxEvaluations, int populationSize, double CR,
				double F, double randomSeed, double **bounds) {
//	clock_t t_ini, t_fin;
	Problem   * problem   ; // The problem to solve
	Algorithm * algorithm ; // The algorithm to use
    Operator  * mutation  ; // Mutation operator

    //problem = new Kursawe("Real", 3);
    //problem = new ZDT6("Real");
//	problem = new DTLZ1("Real", 7, 2);
    problem = new BenchmarkProblem("Real", nreal, nobj, bounds);
	//problem = new DTLZ1("Real");
	//problem = new DTLZ3("Real", 12, 2);
	//problem = new DTLZ4("Real", 12, 2);
    //cout << "El n�mero de objetivos es " << problem->getNumberOfObjectives() << endl;

    //QualityIndicator indicators ; // Object to get quality indicators

    algorithm = new paes(problem);

// Algorithm parameters
	int archiveSizeValue = populationSize;
	int *archiveSizePtr = &archiveSizeValue;
	int biSectionsValue = 5;
	int *biSectionsPtr = &biSectionsValue;
	int maxEvaluationsValue = maxEvaluations;
	int *maxEvaluationsPtr = &maxEvaluationsValue;

   algorithm->setInputParameter("archiveSize",archiveSizePtr);
   algorithm->setInputParameter("biSections",biSectionsPtr);
   algorithm->setInputParameter("maxEvaluations",maxEvaluationsPtr);

   // Mutation (Real variables)
   map<string, void *> parameters; // Operator parameters
   parameters.clear();
   double probabilityValue1 = 1.0/problem->getNumberOfVariables();
   double *probabilityPtr1 = &probabilityValue1;
   double distributionIndexValue1 = 20.0;
   double *distributionIndexPtr1 = &distributionIndexValue1;
   parameters["probability"] =  probabilityPtr1 ;
   parameters["distributionIndex"] = distributionIndexPtr1 ;
   mutation = new PolynomialMutation(parameters);

   // Mutation (BinaryReal variables)
   //mutation = MutationFactory.getMutationOperator("BitFlipMutation");
   //mutation.setParameter("probability",0.1);

   algorithm->addOperator("mutation",mutation);

   // Checking parameters
   //cout << "PolynomialMutation:" << endl;
   //cout << *(double *) algorithm->getOperator("mutation")->getParameter("probability") << endl;
   //cout << *(double *) algorithm->getOperator("mutation")->getParameter("distributionIndex") << endl;


//	t_ini = clock();
	SolutionSet * population = algorithm->execute();
	for (int i = 0; i < population->size(); i++) {
		Solution *solution = population->get(i);
		for (int j = 0; j < nreal; j++) {
			Real *var = (Real *) solution->getDecisionVariables()[j];
			xs[i][j] = var->getValue();
		}
		for (int j = 0; j < nobj; j++)
			fxs[i][j] = solution->getObjective(j);
	}
	
	
	return population->size();
	
	
//	t_fin = clock();
//	double secs = (double) (t_fin - t_ini);
//	secs = secs / CLOCKS_PER_SEC;


	// Result messages
//	cout << "Total execution time: " << secs << "s" << endl;
//	cout << "Variables values have been writen to file VAR" << endl;
//	population->printVariablesToFile("VAR");
//	cout << "Objectives values have been writen to file FUN" << endl;
//	population->printObjectivesToFile("FUN");


 /*
  if (indicators != null) {
	logger_.info("Quality indicators") ;
	logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
	logger_.info("GD         : " + indicators.getGD(population)) ;
	logger_.info("IGD        : " + indicators.getIGD(population)) ;
	logger_.info("Spread     : " + indicators.getSpread(population)) ;
	logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;
  } // if
  */


  // return 10;
} // main

#endif
