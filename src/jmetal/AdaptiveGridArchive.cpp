#include "../config.h"

#ifdef JMETAL

//  AdaptiveGridArchive.cpp
//
//  Author:
//       Cristian Zambrano <cristian_uteq@hotmail.com>
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

#include "AdaptiveGridArchive.h"


AdaptiveGridArchive::AdaptiveGridArchive (int maxSize,int bisections,
		int objectives) : Archive(maxSize) {

	maxSize_   = maxSize;
	dominance_ = new DominanceComparator();
	grid_      = new AdaptiveGrid(bisections,objectives);
}

bool AdaptiveGridArchive::add(Solution * solution, int eval) {
	int i=0;
	Solution * aux; //Store an solution temporally

	while (i < solutionsList_.size()){
	    aux = solutionsList_[i];
	    int flag = dominance_->compare(solution,aux);
	 	if (flag == -1) { // The Individual to insert dominates other
	 		solutionsList_.erase (solutionsList_.begin()+i); //Delete it from the archive
	 		int location = grid_->location(aux);
			if (grid_->getLocationDensity(location) > 1) {//The hypercube contains
				grid_->removeSolution(location);            //more than one individual
			} else {
				grid_->updateGrid(this);
			} // else
		} // if
		else if (flag == 1) { // An Individual into the file dominates the
			// solution to insert
			return false; // The solution will not be inserted
		} // else if
	 	i++;
	} // while

	// At this point, the solution may be inserted
	if (size() == 0){ //The archive is empty
		solutionsList_.push_back(solution);
		grid_->updateGrid(this);
		return true;
	} //

	if (size() < maxSize_){ //The archive is not full
		//cout << "Agregar " << eval << endl;
		grid_->updateGrid(solution,this,eval); // Update the grid if applicable
		//cout << "OK Update " << eval << endl;
		int location ;
		location= grid_->location(solution); // Get the location of the solution
		grid_->addSolution(location); // Increment the density of the hypercube
		//cout << "OK AddSol " << eval << endl;
		solutionsList_.push_back(solution);// Add the solution to the list
		//cout << "OK Agregado " << eval << endl;
		return true;
	} // if

	// At this point, the solution has to be inserted and the archive is full
	grid_->updateGrid(solution,this,eval);
	int location = grid_->location(solution);
	if (location == grid_->getMostPopulated()) { // The solution is in the
		// most populated hypercube
		return false; // Not inserted
	}
	else {
		// Remove an solution from most populated area
		i=0;
		bool removed = false;

		//cout << "Empieza a remover" << endl;
		while (i < solutionsList_.size()){
			if (!removed) {
				aux = solutionsList_[i];
				int location2 = grid_->location(aux);
				if (location2 == grid_->getMostPopulated()) {
					//cout << "Va a remover" << endl;
					solutionsList_.erase (solutionsList_.begin()+i);
					//cout << "Removio " << endl;
					grid_->removeSolution(location2);
				} // if
			} // if
			i++;
		} // while
		// A solution from most populated hypercube has been removed,
		// insert now the solution
		grid_->addSolution(location);
		solutionsList_.push_back(solution);
	} // else
	return true;
} // add

AdaptiveGrid * AdaptiveGridArchive::getGrid() {
	return grid_;
} // AdaptativeGrid



#endif
