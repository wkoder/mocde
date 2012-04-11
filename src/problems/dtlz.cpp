/* Test problem definitions */

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "dtlz.h"
#include "../moead/global.h"

///* Test problem DTLZ1
// * # of real variables = M + k - 1
// * # of binary variables = 0
// * # of objectives = M
// * # of constraints = 0
// */
//void dtlz1(double *xreal, double *obj) {
//	int i = 0;
//	int j = 0;
//	int n = nreal;
//	int k = n - nobj + 1;
//
//	double g = 0;
//	for (i = n - k + 1; i <= n; i++)
//		g += pow(xreal[i - 1] - 0.5, 2) - cos(20 * M_PI * (xreal[i - 1] - 0.5));
//
//	g = 100 * (k + g);
//
//	for (i = 1; i <= nobj; i++) {
//		double f = 0.5 * (1 + g);
//		for (j = nobj - i; j >= 1; j--)
//			f *= xreal[j - 1];
//
//		if (i > 1)
//			f *= 1 - xreal[(nobj - i + 1) - 1];
//
//		obj[i - 1] = f;
//
//		obj[i - 1] += g * 0.00001;
//	}
//}
//
///* Test problem DTLZ2
// * # of real variables = M + k - 1
// * # of binary variables = 0
// * # of objectives = M
// * # of constraints = 0
// */
//void dtlz2(double *xreal, double *obj) {
//	int M = nobj;
//	/*
//	 int K = nreal - M + 1;
//	 int nVars = M + K - 1;
//	 */
//	int nVars = nreal;
//	double g;
//	int i, m; /*, nVars = M + K - 1;*/
//
//	g = 0;
//	for (i = M - 1; i < nVars; ++i)
//		g += pow(xreal[i] - 0.5, 2);
//	//g *= 4;
//
//	for (m = 0; m < M; ++m) {
//		obj[m] = 1.0 + g;
//
//		for (i = 0; i < M - m - 1; ++i)
//			obj[m] *= cos(xreal[i] * M_PI / 2.0);
//
//		if (m > 0)
//			obj[m] *= sin(xreal[M - m - 1] * M_PI / 2.0);
//
////      obj[m] += g*0.01;
//	}
//}
//
//#define a 1
//#define b 10
//#define c 8
//#define c1 (a+c/2)
//#define c2 (c+2*a)
//#define b1 (b/2)
//
//void r_dtlz2(double *x, double *f) {
//	int nx = nreal;
//	int n_obj = nobj;
//	///////////////////////
//
//	int i = 0, j = 0;
//	int k = nx - n_obj + 1;
//	double g = 0;
//	//double z[nx],zz[nx],p[nx],psum[n_obj],M[nx][nx],lamda_l[nx];
//	double *z, *zz, *p, *psum, **M, *lamda_l;
//	double M_10D[10][10] = { { 0.0346, -0.7523, 0.3561, -0.2958, 0.4675, 0, 0, 0, 0, 0 }, { 0.8159, -0.0423, 0.4063, 0.3455, -0.2192, 0, 0,
//			0, 0, 0 }, { -0.3499, 0.3421, 0.8227, -0.2190, -0.1889, 0, 0, 0, 0, 0 }, { -0.0963, -0.4747, -0.0998, -0.2429, -0.8345, 0, 0, 0,
//			0, 0 }, { -0.4487, -0.2998, 0.1460, 0.8283, -0.0363, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 1,
//			0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } };
//	double lamda_l_10D[10] = { 0.313, 0.312, 0.321, 0.316, 0.456, 1, 1, 1, 1, 1 };
//
//	double M_30D[30][30] = { { 0.0128, 0.2165, 0.4374, -0.0800, 0.0886, -0.2015, 0.1071, 0.2886, 0.2354, 0.2785, -0.1748, 0.2147, 0.1649,
//			-0.3043, 0.5316, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 },
//			{ 0.4813, 0.2420, -0.3663, -0.0420, -0.0088, -0.4945, -0.3073, 0.1990, 0.0441, -0.0627, 0.0191, 0.3880, -0.0618, -0.0319,
//					-0.1833, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { 0.4816, -0.2254, 0.0663, 0.4801, 0.2009, -0.0008,
//					-0.1501, 0.0269, -0.2037, 0.4334, -0.2157, -0.3175, -0.0923, 0.1451, 0.1118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0.000, 0.000 }, { -0.0876, -0.2667, -0.0063, 0.2114, 0.4506, 0.0823, -0.0125, 0.2313, 0.0840, -0.2376, 0.1938, -0.0030,
//					0.3391, 0.0863, 0.1231, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { -0.1025, 0.4011, -0.0117, 0.2076,
//					0.2585, 0.1124, -0.0288, 0.3095, -0.6146, -0.2376, 0.1938, -0.0030, 0.3391, 0.0863, 0.1231, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0.000, 0.000 }, { 0.4543, -0.2761, -0.2985, -0.2837, 0.0634, 0.1070, 0.2996, -0.2690, -0.1634, -0.1452,
//					0.1799, -0.0014, 0.2394, -0.2745, 0.3969, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { -0.1422, -0.4364,
//					0.0751, -0.2235, 0.3966, -0.0252, 0.0908, 0.0477, -0.2254, 0.1801, -0.0552, 0.5770, -0.0396, 0.3765, -0.0522, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { 0.3542, -0.2245, 0.3497, -0.1609, -0.1107, 0.0079, 0.2241, 0.4517,
//					0.1309, -0.3355, -0.1123, -0.1831, 0.3000, 0.2045, -0.3191, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, {
//					0.0005, 0.0377, -0.2808, -0.0641, 0.1316, 0.2191, 0.0207, 0.3308, 0.4117, 0.3839, 0.5775, -0.1219, 0.1192, 0.2435,
//					0.0414, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { -0.1177, -0.0001, -0.1992, -0.4533, 0.4234, -0.0191,
//					-0.3740, 0.1325, 0.0972, -0.2042, -0.3493, -0.4018, -0.1087, 0.0918, 0.2217, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0.000, 0.000 }, { 0.1818, 0.3022, -0.1388, -0.2380, -0.0773, 0.6463, 0.0450, 0.1030, -0.0958, 0.2837, -0.3969, 0.1779,
//					-0.0251, -0.1543, -0.2452, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000, 0.000 }, { -0.1889, -0.4397, -0.2206, 0.0981,
//					-0.5203, 0.1325, -0.3427, 0.4242, -0.1271, -0.0291, -0.0795, 0.1213, 0.0565, -0.1092, 0.2720, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0 }, { -0.1808, -0.0624, -0.2689, 0.2289, 0.1128, -0.0844, -0.0549, -0.2202, 0.2450, 0.0825, -0.3319,
//					0.0513, 0.7523, 0.0043, -0.1472, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { -0.0983, 0.0611, -0.4145, 0.3017,
//					0.0410, -0.0703, 0.6250, 0.2449, 0.1307, -0.1714, -0.3045, 0.0218, -0.2837, 0.1408, 0.1633, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0 }, { 0.2026, 0.0324, 0.1496, 0.3129, 0.1437, 0.4331, -0.2629, -0.1498, 0.3746, -0.4366, 0.0163, 0.3316,
//					-0.0697, 0.1833, 0.2412, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
//					0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 1, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 }, {
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//					0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 } };
//
//	double lamda_l_30D[30] = { 0.113, 0.105, 0.117, 0.119, 0.108, 0.110, 0.101, 0.107, 0.111, 0.109, 0.120, 0.108, 0.101, 0.105, 0.116,
//			1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 };
//
//	z = new double[nx];
//	zz = new double[nx];
//	p = new double[nx];
//	psum = new double[n_obj];
//	M = new double*[nx];
//	for (i = 0; i < nx; i++)
//		M[i] = new double[nx];
//	lamda_l = new double[nx];
//
//	if (nx == 10) {
//		for (i = 0; i < nx; i++) {
//			for (j = 0; j < nx; j++) {
//				M[i][j] = M_10D[i][j];
//			}
//			lamda_l[i] = lamda_l_10D[i];
//		}
//	} else {
//		for (i = 0; i < nx; i++) {
//			for (j = 0; j < nx; j++) {
//				M[i][j] = M_30D[i][j];
//			}
//			lamda_l[i] = lamda_l_30D[i];
//		}
//	}
//
//	for (i = 0; i < nx; i++) {
//		z[i] = 0;
//		for (j = 0; j < nx; j++) {
//			z[i] += M[i][j] * x[j];
//		}
//		if (z[i] >= 0 && z[i] <= 1) {
//			zz[i] = z[i];
//			p[i] = 0;
//		} else if (z[i] < 0) {
//			zz[i] = -lamda_l[i] * z[i];
//			p[i] = -z[i];
//		} else {
//			zz[i] = 1 - lamda_l[i] * (z[i] - 1);
//			p[i] = z[i] - 1;
//		}
//	}
//	for (j = 0; j < n_obj; j++) {
//		psum[j] = 0;
//	}
//
//	for (i = nx - k + 1; i <= nx; i++) {
//		g += pow(zz[i - 1] - 0.5, 2);
//		for (j = 0; j < n_obj; j++) {
//			psum[j] = sqrt(pow(psum[j], 2) + pow(p[i - 1], 2));
//		}
//	}
//
//	for (i = 1; i <= n_obj; i++) {
//		double ff = (1 + g);
//		for (j = n_obj - i; j >= 1; j--) {
//			ff *= cos(zz[j - 1] * M_PI / 2.0);
//			psum[i - 1] = sqrt(pow(psum[i - 1], 2) + pow(p[j - 1], 2));
//		}
//		if (i > 1) {
//			ff *= sin(zz[(n_obj - i + 1) - 1] * M_PI / 2.0);
//			psum[i - 1] = sqrt(pow(psum[i - 1], 2) + pow(p[(n_obj - i + 1) - 1], 2));
//		}
//
//		f[i - 1] = 2.0 / (1 + exp(-psum[i - 1])) * (ff + 1);
//	}
//
//	delete[] z;
//	delete[] zz;
//	delete[] p;
//	delete[] psum;
//	delete[] lamda_l;
//	for (i = 0; i < nx; i++)
//		delete[] M[i];
//	delete[] M;
//}
//
///* Test problem DTLZ3
// * # of real variables = M + k - 1
// * # of binary variables = 0
// * # of objectives = M
// * # of constraints = 0
// */
//void dtlz3(double *xreal, double *obj) {
//	int M = nobj;
//	int nVars = nreal;
//	int k = nreal - nobj + 1;
//	double g;
//	int i, m; /*, nVars = M + K - 1;*/
//
//	g = 0;
//	for (i = M - 1; i < nVars; ++i)
//		g += pow(xreal[i] - 0.5, 2) - cos(20 * M_PI * (xreal[i] - 0.5));
//
//	g = 100.0 * (k + g);
//
//	for (m = 0; m < M; ++m) {
//		obj[m] = 1.0 + g;
//
//		for (i = 0; i < M - m - 1; ++i)
//			obj[m] *= cos(xreal[i] * M_PI / 2.0);
//
//		if (m > 0)
//			obj[m] *= sin(xreal[M - m - 1] * M_PI / 2.0);
//
////      obj[m] += g*0.01;
//	}
//}
//
//#define I 3
//
//void dtlz5im(double *x, double *obj) {
//	double g, t;
//	double theta[nobj];
//	int i, m;
//
//	g = 0;
//	for (i = nobj - 1; i < nreal; ++i)
//		g += pow(x[i] - 0.5, 2);
//
//	/* Compute metavariable vector theta */
//	for (i = 0; i < I - 1; i++)
//		theta[i] = x[i] * M_PI / 2.0;
//
//	t = M_PI / (4.0 * (1.0 + g));
//	for (i = I - 1; i < nobj - 1; i++)
//		theta[i] = t * (1 + 2 * g * x[i]);
//
//	for (m = 0; m < nobj; m++) {
//		obj[m] = 1.0 + 100.0 * g;
//
//		for (i = 0; i < nobj - m - 1; i++) {
//			obj[m] *= cos(theta[i]);
//			//if (m == 0)
//			//  printf("t(%d)=%f, cos(t)=%f\n", i, theta[i], cos(theta[i]));
//		}
//
//		if (m > 0)
//			obj[m] *= sin(theta[nobj - m - 1]);
//
//		//obj[m] += 0.0000000001 * g;
//		obj[m] += 15 * g;
//		//obj[m] += 10 * g;
//	}
//}
//
///* Test problem DTLZ7
// * # of real variables = M + k - 1
// * # of binary variables = 0
// * # of objectives = M
// * # of constraints = 0
// */
//void dtlz7(double *xreal, double *obj) {
//	int i;
//	double h, g;
//	int M = nobj;
//	int nVars = nreal;
//
//	/* Compute g */
//	g = 0.0;
//	for (i = M - 1; i < nVars; ++i)
//		g += xreal[i];
//
//	g = 1 + (9.0 / (nVars - M + 1)) * g;
//
//	for (i = 0; i < M - 1; ++i) {
//		obj[i] = xreal[i];
//	}
//
//	/* Compute h */
//	h = 0.0;
//	for (i = 0; i < M - 1; ++i)
//		h += (obj[i] / (1 + g)) * (1 + sin(3 * M_PI * obj[i]));
//	h = M - h;
//
//	obj[M - 1] = (1 + g) * h;
//
//	for (i = 0; i < M - 1; ++i) {
//		obj[i] += g * 0.1;
//	}
//
//	//                    3              4         5           6           7
//	double dtlz7lower[] = { 2.61400874, 2.92101312, 3.22801749, 3.53502186, 3.84202623,
//	//                    8              9        10          11       	  12
//			4.14903061, 4.45603498, 4.76303935, 5.07004372, 5.37704810,
//			//                   13     		14        15
//			5.68405247, 5.99105684, 6.29806121 };
//
//	obj[M - 1] = (obj[M - 1] - dtlz7lower[nobj - 3]) / (2 * nobj - dtlz7lower[nobj - 3]);
//}

void dtlz1(double *xreal, double *obj) {
	int j;
	double sum;
	double g;
	int n = nreal;
	int m = nobj;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++) {
		sum += pow(xreal[j] - 0.5, 2.0) - cos(20.0 * M_PI * (xreal[j] - 0.5));
	}
	g = 100.0 * (k + sum);

	obj[0] = (0.5 * (1.0 + g) * xreal[0] * xreal[1]);
	obj[1] = 0.5 * (1.0 + g) * (1.0 - xreal[1]) * xreal[0];
	obj[2] = 0.5 * (1.0 + g) * (1.0 - xreal[0]);
	return;
}

/* DTLZ2 - MOP */
void dtlz2(double *xreal, double *obj) {
	int j;
	double sum;
	double g;
	int n = nreal;
	int m = nobj;

	sum = 0.0;
	for (j = m - 1; j < n; j++) {
		sum += pow(xreal[j] - 0.5, 2.0);
	}
	g = sum;

	obj[0] = (1.0 + g) * cos(xreal[0] * M_PI * 0.5) * cos(xreal[1] * M_PI * 0.5);
	obj[1] = (1.0 + g) * cos(xreal[0] * M_PI * 0.5) * sin(xreal[1] * M_PI * 0.5);
	obj[2] = (1.0 + g) * sin(xreal[0] * M_PI * 0.5);
	return;
}

/* DTLZ3 - MOP */
void dtlz3(double *xreal, double *obj) {
	int j;
	double sum;
	double g;
	int n = nreal;
	int m = nobj;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++) {
		sum += pow(xreal[j] - 0.5, 2.0) - cos(20.0 * M_PI * (xreal[j] - 0.5));
	}
	g = 100.0 * (k + sum);

	obj[0] = (1.0 + g) * cos(xreal[0] * M_PI * 0.5) * cos(xreal[1] * M_PI * 0.5);
	obj[1] = (1.0 + g) * cos(xreal[0] * M_PI * 0.5) * sin(xreal[1] * M_PI * 0.5);
	obj[2] = (1.0 + g) * sin(xreal[0] * M_PI * 0.5);
	return;
}

/* DTLZ4 - MOP */
void dtlz4(double *xreal, double *obj) {
	int j;
	double alpha;
	double g, factor1, factor2;
	int n = nreal;
	int m = nobj;

	alpha = 100;
	g = 0.0;
	for (j = m - 1; j < n; j++) {
		g += pow(xreal[j] - 0.5, 2.0);
	}

	factor1 = pow(xreal[0], alpha);
	factor2 = pow(xreal[1], alpha);
	obj[0] = (1.0 + g) * cos(factor1 * M_PI * 0.5) * cos(factor2 * M_PI * 0.5);
	obj[1] = (1.0 + g) * cos(factor1 * M_PI * 0.5) * sin(factor2 * M_PI * 0.5);
	obj[2] = (1.0 + g) * sin(factor1 * M_PI * 0.5);
	return;
}

/* DTLZ5 - MOP */
void dtlz5(double *xreal, double *obj) {
	int j;
	double g;
	double tetha[3];
	int n = nreal;
	int m = nobj;

	g = 0.0;
	for (j = m - 1; j < n; j++) {
		g += pow(xreal[j] - 0.5, 2.0);
	}

	tetha[0] = (M_PI * 0.5) * xreal[0];
	tetha[1] = (M_PI / (4.0 * (1.0 + g))) * (1.0 + 2.0 * g * xreal[1]);

	obj[0] = (1.0 + g) * cos(tetha[0]) * cos(tetha[1]);
	obj[1] = (1.0 + g) * cos(tetha[0]) * sin(tetha[1]);
	obj[2] = (1.0 + g) * sin(tetha[0]);
	return;
}

/* DTLZ6 - MOP */
void dtlz6(double *xreal, double *obj) {
	int j;
	double g;
	double tetha[3];
	int n = nreal;
	int m = nobj;

	g = 0.0;
	for (j = m - 1; j < n; j++) {
		g += pow(xreal[j], 0.1);
	}

	tetha[0] = (M_PI * 0.5) * xreal[0];
	tetha[1] = (M_PI / (4.0 * (1.0 + g))) * (1.0 + 2.0 * g * xreal[1]);
	obj[0] = (1.0 + g) * cos(tetha[0]) * cos(tetha[1]);
	obj[1] = (1.0 + g) * cos(tetha[0]) * sin(tetha[1]);
	obj[2] = (1.0 + g) * sin(tetha[0]);
	return;
}

/* DTLZ7 - MOP */
void dtlz7(double *xreal, double *obj) {
	int j;
	double sum;
	double g, h;
	int n = nreal;
	int m = nobj;
	int k = n - m + 1;

	sum = 0.0;
	for (j = m - 1; j < n; j++) {
		sum += xreal[j];
	}
	g = 1.0 + (9.0 / (double) k) * sum;

	obj[0] = xreal[0];
	obj[1] = xreal[1];
	h = m
			- (((obj[0] / (1.0 + g)) * (1.0 + sin(3.0 * M_PI * obj[0])))
					+ ((obj[1] / (1.0 + g)) * (1.0 + sin(3.0 * M_PI * obj[1]))));
	obj[2] = (1.0 + g) * h;
	return;
}
