/* Test problem definitions */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "wfg.h"
#include "../moead/global.h"

using std::vector;

using namespace std;

#define k_factor 1
#define l_factor 10

void wfg1(double *xreal, double *obj) {
//	unsigned i;
	int K = k_factor * (nobj - 1); //# of position-related variables
	int L = l_factor * 2; //# of distance-related variables

	if (nreal != K + L) {
		printf("\nEl numero de variables debe ser %d, y hay %d", K + L, nreal);
		exit(0);
	}

	vector<double> z(xreal, xreal + nreal);
	vector<double> S(nobj);
	for (int m = 0; m < nobj; ++m)
		S[m] = 2 * (m + 1);

//	vector<double> fx = Problems::WFG1(z, K, nobj, S);
//
//	for (i = 0; i < fx.size(); ++i)
//		obj[i] = fx[i];
}

void wfg2(double *xreal, double *obj) {
//	unsigned i;

	int K = k_factor * (nobj - 1); //# of position-related variables
	int L = l_factor * 2; //# of distance-related variables

	if (nreal != K + L) {
		printf("\nEl numero de variables debe ser %d, y hay %d\n\n", (K + L), nreal);
		exit(0);
	}

	vector<double> z(xreal, xreal + nreal);
	vector<double> S(nobj);
	for (int m = 0; m < nobj; ++m)
		S[m] = 1; //2*(m+1);

//	vector<double> fx = Problems::WFG2(z, K, nobj, S);
//
//	for (i = 0; i < fx.size(); ++i)
//		obj[i] = fx[i];
}

void wfg6(double *xreal, double *obj) {
//	unsigned i;
	int K = k_factor * (nobj - 1); //# of position-related variables
	int L = l_factor * 2; //# of distance-related variables

	if (nreal != K + L) {
		printf("\nEl numero de variables debe ser %d, y hay %d\n\n", (K + L), nreal);
		exit(0);
	}

	vector<double> z(xreal, xreal + nreal);
	vector<double> S(nobj);
	for (int m = 0; m < nobj; ++m)
		S[m] = 1; //2*(m+1);

//	vector<double> fx = Problems::WFG6(z, K, nobj, S);
//
//	for (i = 0; i < fx.size(); ++i)
//		obj[i] = fx[i];
}
