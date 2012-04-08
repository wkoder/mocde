#include "../config.h"

#ifdef JMETAL

//  DTLZ1.cpp
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

#include "BenchmarkProblem.h"
#include "../benchmark.h"

BenchmarkProblem::BenchmarkProblem(string solutionType, int numberOfVariables, int numberOfObjectives, double (*bounds)[2]) {
	numberOfVariables_ = numberOfVariables;
	numberOfObjectives_ = numberOfObjectives;
	numberOfConstraints_ = 0;
	problemName_ = "DTLZ1";

	lowerLimit_ = new double[numberOfVariables_]; //(double *)malloc(sizeof(double)*numberOfVariables);
	if (lowerLimit_ == NULL) {
		cout << "Impossible to reserve memory for storing the variable lower limits" << endl;
		exit(-1);
	}

	upperLimit_ = new double[numberOfVariables_]; //(double *)malloc(sizeof(double)*numberOfVariables);
	if (upperLimit_ == NULL) {
		cout << "Impossible to reserve memory for storing the variable lower limits" << endl;
		exit(-1);
	}

	for (int i = 0; i < numberOfVariables_; i++) {
//		lowerLimit_[i] = 0.0;
//		upperLimit_[i] = 1.0;
		lowerLimit_[i] = bounds[i][0];
		upperLimit_[i] = bounds[i][1];
	}

	if (solutionType.compare("BinaryReal") == 0)
		solutionType_ = new BinaryRealSolutionType(this);
	else if (solutionType.compare("Real") == 0) {
		solutionType_ = new RealSolutionType(this);
		cout << "Tipo seleccionado Real" << endl;
	} else if (solutionType.compare("ArrayReal") == 0)
		solutionType_ = new ArrayRealSolutionType(this);
	else {
		cout << "Error: solution type " << solutionType << " invalid" << endl;
		exit(-1);
	}

	fx_ = new double[numberOfObjectives_];
	x_ = new double[numberOfVariables_];
}

BenchmarkProblem::~BenchmarkProblem() {
	delete[] lowerLimit_;
	delete[] upperLimit_;
	delete solutionType_;
	delete[] fx_;
	delete[] x_;
}

/**
 * Evaluates a solution
 * @param solution The solution to evaluate
 */
void BenchmarkProblem::evaluate(Solution *solution) {
	XReal * vars = new XReal(solution);

	for (int i = 0; i < numberOfVariables_; i++)
		x_[i] = vars->getValue(i);

	benchmark::evaluate(x_, fx_);

	for (int i = 0; i < numberOfObjectives_; i++)
		solution->setObjective(i, fx_[i]);

	delete vars;
} // evaluate

#endif
