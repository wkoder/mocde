#ifndef __COMMON_H_
#define __COMMON_H_

#include "global.h"

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}


double dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

bool   dominate(vector<double> &u, vector<double> &v, double epsilon)
{
    int dim = u.size();
	for(int i=0; i<dim; i++)
	{
	    if(u[i]<v[i]-epsilon)
		{
		    return false;
		}
	}
	return true;
}

double fitnessfunction(vector <double> &y_obj, vector <double> &namda)
{
    // Chebycheff Scalarizing Function
	double max_fun = -1.0e+30;

	for(int n=0; n<nobj; n++)
	{
		double diff = fabs(y_obj[n] - idealpoint[n]);
		double feval;
		if(namda[n]==0) 
			feval = 0.0001*diff;
		else
			feval = diff*namda[n];
		if(feval>max_fun) max_fun = feval;

	}

	return max_fun;
}

void load_data(char filename[1024], vector< vector<double> > &data, int dim)
{
	std::ifstream readf(filename);
	vector<double> vect = vector<double>(dim, 0);
	while(!readf.eof())
	{
        for(int i=0; i<dim; i++)
		{
			readf>>vect[i];
			//printf(" %f ", vect[i]);
		}
		data.push_back(vect);
		//vect.clear();    // bug here. 
	}

	readf.close();    
}

void permutation(vector<int> &perm)
{
	unsigned int i, j, size = perm.size();
	vector<int> index(size);
	for(i=0; i<size; i++) index[i] = i;
	for(i=0; i<size; i++)
	{
		j = int(size*rnd_uni(&rnd_uni_init));
		while(1)
		{
			if(index[j]>=0)
			{
				perm[i] = index[j];
				index[j]= -1;
				break;
			}
			else if(j==size-1)
			{
				j=0;
			}
			else
			{
				j++;
			}
		}
	}
}


#endif
