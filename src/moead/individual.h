#include "global.h"

#ifdef MOEAD_IMPL

#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "../benchmark.h"

class CIndividual {
public:
	CIndividual();
	virtual ~CIndividual();

	vector<double> x_var;
	vector<double> y_obj;
	int rank;
	int count;

	void rnd_init();
	void obj_eval();
	void show_objective();
	void show_variable();

	bool operator<(const CIndividual &ind2);
	bool operator<<(const CIndividual &ind2);
	bool operator==(const CIndividual &ind2);
	void operator=(const CIndividual &ind2);
};

CIndividual::CIndividual() {
	x_var = vector<double>(nreal, 0);
	y_obj = vector<double>(nobj, 0);
	rank = 0;
}

CIndividual::~CIndividual() {

}

void CIndividual::rnd_init() {
	double **bounds = benchmark::getBounds();
	for (int n = 0; n < nreal; n++)
		x_var[n] = bounds[n][0] + rnd_uni(&rnd_uni_init) * (bounds[n][1] - bounds[n][0]);
}

void CIndividual::obj_eval() {
	double x[x_var.size()];
	double y[y_obj.size()];
	for (unsigned int i = 0; i < x_var.size(); i++)
		x[i] = x_var[i];
	for (unsigned int i = 0; i < y_obj.size(); i++)
		y[i] = y_obj[i];
	benchmark::evaluate(x, y);
	for (unsigned int i = 0; i < x_var.size(); i++)
		x_var[i] = x[i];
	for (unsigned int i = 0; i < y_obj.size(); i++)
		y_obj[i] = y[i];
//	if(!strcmp("UF1", strTestInstance))  CEC09_F1(y_obj, x_var);
//	if(!strcmp("UF2", strTestInstance))  CEC09_F2(y_obj, x_var);
//	if(!strcmp("UF3", strTestInstance))  CEC09_F3(y_obj, x_var);
//	if(!strcmp("UF4", strTestInstance))  CEC09_F4(y_obj, x_var);
//	if(!strcmp("UF5", strTestInstance))  CEC09_F5(y_obj, x_var);
//	if(!strcmp("UF6", strTestInstance))  CEC09_F6(y_obj, x_var);
//	if(!strcmp("UF7", strTestInstance))  CEC09_F7(y_obj, x_var);
//	if(!strcmp("UF8", strTestInstance))  CEC09_F8(y_obj, x_var);
//	if(!strcmp("UF9", strTestInstance))  CEC09_F9(y_obj, x_var);
//	if(!strcmp("UF10", strTestInstance)) CEC09_F10(y_obj, x_var);
//
//
//	if(!strcmp("R2_DTLZ2_M5", strTestInstance))	CEC09_R2_DTLZ2_M5(y_obj, x_var);	
//	if(!strcmp("R3_DTLZ3_M5", strTestInstance)) CEC09_R3_DTLZ3_M5(y_obj, x_var);
//	if(!strcmp("WFG1_M5", strTestInstance))     CEC09_WFG1_M5(y_obj, x_var);
}

void CIndividual::show_objective() {
	for (int n = 0; n < nobj; n++)
		printf("%f ", y_obj[n]);
	printf("\n");
}

void CIndividual::show_variable() {
	for (int n = 0; n < nreal; n++)
		printf("%f ", x_var[n]);
	printf("\n");
}

void CIndividual::operator=(const CIndividual &ind2) {
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	rank = ind2.rank;
}

bool CIndividual::operator<(const CIndividual &ind2) {
	bool dominated = true;
	for (int n = 0; n < nobj; n++) {
		if (ind2.y_obj[n] < y_obj[n])
			return false;
	}
	if (ind2.y_obj == y_obj)
		return false;
	return dominated;
}

bool CIndividual::operator<<(const CIndividual &ind2) {
	bool dominated = true;
	for (int n = 0; n < nobj; n++) {
		if (ind2.y_obj[n] < y_obj[n] - 0.0001)
			return false;
	}
	if (ind2.y_obj == y_obj)
		return false;
	return dominated;
}

bool CIndividual::operator==(const CIndividual &ind2) {
	if (ind2.y_obj == y_obj)
		return true;
	else
		return false;
}

// defining subproblem 

class CSubproblem {
public:
	CSubproblem();
	virtual ~CSubproblem();

	void show();

	CIndividual indiv; // best solution
	CIndividual saved; // last solution
	vector<double> namda; // weight vector
	vector<int> table; // neighbourhood table

	double fitness;

	void operator=(const CSubproblem &sub2);
};

CSubproblem::CSubproblem() {
	namda = vector<double>(nobj, 0);
}

CSubproblem::~CSubproblem() {
}

void CSubproblem::show() {
	for (unsigned int n = 0; n < namda.size(); n++) {
		printf("%f ", namda[n]);
	}
	printf("\n");
}

void CSubproblem::operator=(const CSubproblem &sub2) {
	indiv = sub2.indiv;
	table = sub2.table;
	namda = sub2.namda;
}

#endif

#endif
