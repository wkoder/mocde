#include "../config.h"

#ifdef JMETAL

//  OverallConstraintViolationComparator.cpp
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


#include "OverallConstraintViolationComparator.h"


/**
 * Compares two solutions.
 * @param o1 Object representing the first <code>Solution</code>.
 * @param o2 Object representing the second <code>Solution</code>.
 * @return -1, or 0, or 1 if o1 is less than, equal, or greater than o2,
 * respectively.
 */
int OverallConstraintViolationComparator::compare(Solution *o1, Solution *o2) {
  double overall1, overall2;
  overall1 = o1->getOverallConstraintViolation();
  overall2 = o2->getOverallConstraintViolation();

  if ((overall1 < 0) && (overall2 < 0)) {
    if (overall1 > overall2){
      return -1;
    } else if (overall2 > overall1){
      return 1;
    } else {
      return 0;
    }
  } else if ((overall1 == 0) && (overall2 < 0)) {
    return -1;
  } else if ((overall1 < 0) && (overall2 == 0)) {
    return 1;
  } else {
    return 0;
  }
} // compare


#endif
