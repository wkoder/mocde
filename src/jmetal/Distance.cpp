#include "../config.h"

#ifdef JMETAL

//  Distance.cpp
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
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


#include "Distance.h"


/**
 * This class implements some utilities for calculating distances
 */


/**
 * Constructor.
 */
Distance::Distance() {
  //do nothing.
} // Distance


/**
 * Returns a matrix with distances between solutions in a
 * <code>SolutionSet</code>.
 * @param solutionSet The <code>SolutionSet</code>.
 * @return a matrix with distances.
 */
double ** Distance::distanceMatrix(SolutionSet * solutionSet) {
  Solution * solutionI;
  Solution * solutionJ;

  //The matrix of distances
  int size = solutionSet->size();
  double ** distance;// = new double [size][size];
  for (int i = 0; i < solutionSet->size(); i++){
    distance[i] = new double[size];
  }
  //-> Calculate the distances
  for (int i = 0; i < solutionSet->size(); i++){
    distance[i][i] = 0.0;
    solutionI = solutionSet->get(i);
    for (int j = i + 1; j < solutionSet->size(); j++){
      solutionJ = solutionSet->get(j);
      distance[i][j] = this->distanceBetweenObjectives(solutionI,solutionJ);
      distance[j][i] = distance[i][j];
    } // for
  } // for

  //->Return the matrix of distances
  return distance;
} // distanceMatrix


/**
 * Returns the minimum distance from a <code>Solution</code> to a
 * <code>SolutionSet according to the objective values</code>.
 * @param solution The <code>Solution</code>.
 * @param solutionSet The <code>SolutionSet</code>.
 * @return The minimum distance between solution and the set.
 */
double Distance::distanceToSolutionSetInObjectiveSpace(Solution *solution,
                                    SolutionSet * solutionSet) {
  //At start point the distance is the max
  double distance = std::numeric_limits<double>::max();

  // found the min distance respect to population
  for (int i = 0; i < solutionSet->size();i++){
    double aux = this->distanceBetweenObjectives(solution,solutionSet->get(i));
    if (aux < distance)
      distance = aux;
  } // for

  //->Return the best distance
  return distance;
} // distanceToSolutionSetinObjectiveSpace


/**
 * Returns the minimum distance from a <code>Solution</code> to a
 * <code>SolutionSet according to the variable values</code>.
 * @param solution The <code>Solution</code>.
 * @param solutionSet The <code>SolutionSet</code>.
 * @return The minimum distance between solution and the set.
 */
double Distance::distanceToSolutionSetInSolutionSpace(Solution * solution,
                                  SolutionSet * solutionSet) {
  //At start point the distance is the max
  double distance = std::numeric_limits<double>::max();

  // found the min distance respect to population
  for (int i = 0; i < solutionSet->size(); i++) {
    double aux = this->distanceBetweenSolutions(solution,solutionSet->get(i));
    if (aux < distance)
      distance = aux;
  } // for

  //->Return the best distance
  return distance;
} // distanceToSolutionSetInSolutionSpace


/**
 * Returns the distance between two solutions in the search space.
 * @param solutionI The first <code>Solution</code>.
 * @param solutionJ The second <code>Solution</code>.
 * @return the distance between solutions.
 */
double Distance::distanceBetweenSolutions(Solution * solutionI, Solution * solutionJ) {
  //->Obtain his decision variables
  Variable ** decisionVariableI = solutionI->getDecisionVariables();
  Variable ** decisionVariableJ = solutionJ->getDecisionVariables();

  double diff;    //Auxiliar var
  double distance = 0.0;
  //-> Calculate the Euclidean distance
  for (int i = 0; i < solutionI->getNumberOfVariables(); i++){
    diff = decisionVariableI[i]->getValue() -
           decisionVariableJ[i]->getValue();
    distance += pow(diff,2.0);
  } // for

  //-> Return the euclidean distance
  return sqrt(distance);
} // distanceBetweenSolutions


/**
 * Returns the distance between two solutions in objective space.
 * @param solutionI The first <code>Solution</code>.
 * @param solutionJ The second <code>Solution</code>.
 * @return the distance between solutions in objective space.
 */
double Distance::distanceBetweenObjectives(Solution * solutionI, Solution * solutionJ) {
  double diff;    //Auxiliar var
  double distance = 0.0;
  //-> Calculate the euclidean distance
  for (int nObj = 0; nObj < solutionI->getNumberOfObjectives(); nObj++){
    diff = solutionI->getObjective(nObj) - solutionJ->getObjective(nObj);
    distance += pow(diff,2.0);
  } // for

  //Return the euclidean distance
  return sqrt(distance);
} // distanceBetweenObjectives


/**
 * Assigns crowding distances to all solutions in a <code>SolutionSet</code>.
 * @param solutionSet The <code>SolutionSet</code>.
 * @param nObjs Number of objectives.
 */
void Distance::crowdingDistanceAssignment(SolutionSet * solutionSet, int nObjs) {
  int size = solutionSet->size();

  if (size == 0)
    return;

  if (size == 1) {
    solutionSet->get(0)->setCrowdingDistance(std::numeric_limits<double>::max());
    return;
  } // if

  if (size == 2) {
    solutionSet->get(0)->setCrowdingDistance(std::numeric_limits<double>::max());
    solutionSet->get(1)->setCrowdingDistance(std::numeric_limits<double>::max());
    return;
  } // if

  //Use a new SolutionSet to evite alter original solutionSet
  SolutionSet * front = new SolutionSet(size);
  for (int i = 0; i < size; i++){
    front->add(solutionSet->get(i));
  }

  for (int i = 0; i < size; i++)
    front->get(i)->setCrowdingDistance(0.0);

  double objetiveMaxn;
  double objetiveMinn;
  double distance;

  for (int i = 0; i<nObjs; i++) {
    // Sort the population by Obj n
    Comparator * c = new ObjectiveComparator(i);
    front->sort(c);
    delete c;
    objetiveMinn = front->get(0)->getObjective(i);
    objetiveMaxn = front->get(front->size()-1)->getObjective(i);

    //Set de crowding distance
    front->get(0)->setCrowdingDistance(std::numeric_limits<double>::max());
    front->get(size-1)->setCrowdingDistance(std::numeric_limits<double>::max());

    for (int j = 1; j < size-1; j++) {
      distance = front->get(j+1)->getObjective(i) - front->get(j-1)->getObjective(i);
      distance = distance / (objetiveMaxn - objetiveMinn);
      distance += front->get(j)->getCrowdingDistance();
      front->get(j)->setCrowdingDistance(distance);
    } // for
  } // for

  front->clear();
  delete front;

} // crowdingDistanceAssignment


#endif
