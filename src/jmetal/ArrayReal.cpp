#include "../config.h"

#ifdef JMETAL

//  ArrayReal.cpp
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


#include "ArrayReal.h"


/**
 * Constructor
 */
ArrayReal::ArrayReal() {
  problem_ = NULL;
  size_   = 0;
  array_ = NULL;
} // Constructor


/**
 * Constructor
 * @param size Size of the array
 */
ArrayReal::ArrayReal(int size, Problem * problem) {
  problem_ = problem;
  size_   = size;
  array_ = new double[size_];

  for (int i = 0; i < size_ ; i++) {
    array_[i] = PseudoRandom::randDouble()*(problem_->getUpperLimit(i)-
                                           problem_->getLowerLimit(i))+
                                           problem_->getLowerLimit(i);
  } // for
} // Constructor


/**
 * Copy Constructor
 * @param arrayReal The arrayDouble to copy
 */
ArrayReal::ArrayReal(ArrayReal * arrayReal) {
  problem_ = arrayReal->problem_ ;
  size_  = arrayReal->size_;
  array_ = new double[size_];

  for (int i = 0; i < size_; i++) {
    array_[i] = arrayReal->array_[i];
  } // for
} // Copy Constructor


/**
 * Destructor
 */
ArrayReal::~ArrayReal() {
  delete [] array_;
} // ~ArrayReal


/**
 * Creates an exact copy of a <code>BinaryReal</code> object.
 * @return The copy of the object
 */
Variable * ArrayReal::deepCopy() {
  return new ArrayReal(this);
} // deepCopy


/**
 * Returns the length of the arrayReal.
 * @return The length
 */
int ArrayReal::getLength(){
  return size_;
} // getLength


/**
 * getValue
 * @param index Index of value to be returned
 * @return the value in position index
 */
double ArrayReal::getValue(int index) {
  if ((index >= 0) && (index < size_))
    return array_[index] ;
  else {
    cout << "ArrayReal.getValue: index value (" << index << ") invalid" << endl;
  } // if
} // getValue

double ArrayReal::getValue() {
  cout << "ERROR: ArrayReal::getValue() without index" << endl;
} // getValue


/**
 * setValue
 * @param index Index of value to be returned
 * @param value The value to be set in position index
 */
void ArrayReal::setValue(int index, double value) {
  if ((index >= 0) && (index < size_))
    array_[index] = value;
  else {
    cout << "ArrayReal.getValue: index value (" << index << ") invalid" << endl;
  } // if
} // setValue

void ArrayReal::setValue(double value) {
  cout << "ERROR: ArrayReal::setValue(value) without index" << endl;
} // setValue


/**
 * Get the lower bound of a value
 * @param index The index of the value
 * @return the lower bound
 */
double ArrayReal::getLowerBound(int index) {
  if ((index >= 0) && (index < size_))
    return problem_->getLowerLimit(index) ;
  else {
    cout << "ArrayReal.getValue: index value (" << index << ") invalid" << endl;
  } // if
} // getLowerBound

double ArrayReal::getLowerBound() {
  cout << "ERROR: ArrayReal::getLowerBound() without index" << endl;
} // getLowerBound


/**
 * Get the upper bound of a value
 * @param index The index of the value
 * @return the upper bound
 */
double ArrayReal::getUpperBound(int index) {
  if ((index >= 0) && (index < size_))
    return problem_->getUpperLimit(index);
  else {
    cout << "ArrayReal.getValue: index value (" << index << ") invalid" << endl;
  } // if
} // getLowerBound

double ArrayReal::getUpperBound() {
  cout << "ERROR: ArrayReal::getUpperBound() without index" << endl;
} // getUpperBound


/**
 * Returns a string representing the object
 * @return The string
 */
string ArrayReal::toString(){
  string string_ = "";
  for (int i = 0; i < (size_ - 1); i ++) {
    string_ += array_[i];
    string_ += " ";
  }
  string_ += array_[size_ -1];
  return string_ ;
} // toString

#endif
