/*==========================================================================
//  C++ Implementation of MOEA/D Based on Differential Evolution (DE) for Contest Multiobjective
//  Problems in CEC2009
//
//  Author: Hui Li
//
//  See the details of MOEA/D-DE and test problems in the following papers
//
//  1) H. Li and Q. Zhang, Comparison Between NSGA-II and MOEA/D on a Set of Multiobjective Optimization
//  Problems with Complicated Pareto Sets, Technical Report CES-476, Department of Computer Science,
//  University of Essex, 2007
//
//  2) H. Li and Q. Zhang, Multiobjective Optimization Problems with Complicated Pareto Sets, MOEA/D and NSGA-II,
//  IEEE Transaction on Evolutionary Computation, 2008, to appear.
//
//  If you have any questions about the codes, please contact
//  Dr. Hui Li       at   hzl@cs.nott.ac.uk   or
//  Dr. Qingfu Zhang at   qzhang@essex.ac.uk
//
//  Date: 14/NOV/2008
//
// ===========================================================================*/
#include "time.h"
#include "global.h"
#include "algorithm.h"

int main(int argc,char *argv[])
{
	int ss;
	for(ss=1; ss<argc; ss++)
	{
		// random number init
		rnd_uni_init = 90.0;

		std::ifstream readf("TestInstance.txt");

		int numOfInstance;
		readf>>numOfInstance;

		printf("-- % Instances are being tested with seed %d---\n", atoi(argv[ss]));

		for(int inst=1; inst<=13; inst++)
		{
			// the parameter setting of test instance
			readf>>strTestInstance;
			readf>>nreal;
			readf>>nobj;

			printf("-- Instance: %s, %d variables %d objectives \n", strTestInstance, nreal, nobj);

			for(int run=1; run<=30; run++)
			{
				printf("%d ", run);
				CMOEAD MOEAD;
				MOEAD.load_parameter();
				MOEAD.exec_emo(run);
			}
			printf(endl);
		}
	}
	return 0;
}
