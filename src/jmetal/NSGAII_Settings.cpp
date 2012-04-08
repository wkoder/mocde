#include "../config.h"

#ifdef JMETAL

//  NSGAII_Settings.cpp
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
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

#include "NSGAII_Settings.h"

/**
 * Default constructor
 */
NSGAII_Settings::NSGAII_Settings () {
} // Settings

/**
 * Constructor
 */
NSGAII_Settings::NSGAII_Settings(char * problemName) {
	problemName_ = problemName ;

  problem_ = ProblemFactory::getProblem(problemName_);

  // Algorithm parameters
  populationSize_ = 100;
  maxEvaluations_ = 25000;
  mutationProbability_         = 1.0/problem_->getNumberOfVariables() ;
  crossoverProbability_        = 0.9   ;
  mutationDistributionIndex_   = 20.0  ;
  crossoverDistributionIndex_  = 20.0  ;
} // Settings

/**
 * Configure method
 */
Algorithm * NSGAII_Settings::configure() {
	Algorithm * algorithm ;
  Operator  * crossover ; // Crossover operator
  Operator  * mutation  ; // Mutation operator
  Operator  * selection ; // Selection operator

	algorithm = new NSGAII(problem_);
  algorithm->setInputParameter("populationSize",&populationSize_);
  algorithm->setInputParameter("maxEvaluations",&maxEvaluations_);

	// Mutation and Crossover for Real codification
	map<string, void *> parameters;

  double crossoverProbability = crossoverProbability_;
  double crossoverDistributionIndex = crossoverDistributionIndex_ ;
  parameters["probability"] =  &crossoverProbability;
  parameters["distributionIndex"] = &crossoverDistributionIndex;
  crossover = new SBXCrossover(parameters);

	parameters.clear();
  double mutationProbability = mutationProbability_;
  double mutationDistributionIndex = mutationDistributionIndex_;
  parameters["probability"] = &mutationProbability;
  parameters["distributionIndex"] = &mutationDistributionIndex;
  mutation = new PolynomialMutation(parameters);

	// Selection Operator
	parameters.clear();
	selection = new BinaryTournament2(parameters);

	// Add the operators to the algorithm
	algorithm->addOperator("crossover",crossover);
	algorithm->addOperator("mutation",mutation);
	algorithm->addOperator("selection",selection);

	cout << "NGSAII algorithm initialized." << endl;

	return algorithm ;
}


#endif
