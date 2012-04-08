#include "../config.h"

#ifdef JMETAL

//  RealSolutionType.h
//
//  Author:
//       Juan J. Durillo <durillo@lcc.uma.es>
//       Esteban L�pez <esteban@lcc.uma.es>
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

#ifndef __REAL_SOLUTION_TYPE__
#define __REAL_SOLUTION_TYPE__

#include "SolutionType.h"
#include "Real.h"
#include <string>

/**
 * This class is aimed at defining a Type encoding a Real solution
**/
class RealSolutionType : public SolutionType {

public:
  RealSolutionType(Problem *problem);
  Variable **createVariables();

};

#endif

#endif
