#include "../config.h"

#ifdef JMETAL

#ifndef _BENCHMARK_PROBLEM_H_
#define _BENCHMARK_PROBLEM_H_

#include "Problem.h"
#include <math.h>
#include "BinaryRealSolutionType.h"
#include "RealSolutionType.h"
#include "ArrayRealSolutionType.h"
#include "XReal.h"
#include "Solution.h"

class BenchmarkProblem : public Problem {
public:
	BenchmarkProblem(string solutionType, int numberOfVariables = 7, int numberOfObjectives = 3, double **bounds = NULL);
	void evaluate(Solution *solution);

	virtual ~BenchmarkProblem();
private:
	double * fx_ ;
	double * x_  ;
};


#endif
#endif
