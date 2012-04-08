#include "../config.h"

#ifdef JMETAL

//  ZDT6.h
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

#ifndef __ZDT6_H__
#define __ZDT6_H__

#include "Problem.h"
#include <math.h>
#include <cmath>
#include "BinaryRealSolutionType.h"
#include "RealSolutionType.h"
#include "ArrayRealSolutionType.h"
#include "XReal.h"
#include "Solution.h"

class ZDT6 : public Problem {
private:
	double evalG(XReal *x) ;
	double evalH(double f, double g) ;
	double * fx_ ;

public:
	ZDT6(string solutionType, int numberOfVariables = 10);
	void evaluate(Solution *solution);

	virtual ~ZDT6();
};

#endif /* __ZDT6_H__ */

#endif
