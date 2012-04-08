#include "../config.h"

#ifdef JMETAL

//  ssGA.h
//
//  Author:
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

#ifndef __SSGA_H__
#define __SSGA_H__

#include "Algorithm.h"
#include "Problem.h"
#include "SolutionSet.h"
#include "ObjectiveComparator.h"
#include "WorstSolutionSelection.h"

/**
 * Class implementing a steady-state genetic algorithm
 */
class ssGA : public Algorithm {

public:
  ssGA(Problem * problem);
  SolutionSet * execute();

};

#endif /* __SSGA_H__ */

#endif
