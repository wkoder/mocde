/*
 * stats.h
 *
 *  Created on: Jun 15, 2012
 *      Author: Moises Osorio
 */

#ifndef STATS_H_
#define STATS_H_

#include <stdlib.h>
#include <string>
#include <vector>

#include "individual.h"

namespace stats {
	void configure(int outputInterval, std::string prefix);
	bool timeToReport();
	void report(std::vector<Individual *> population);
}

#endif /* STATS_H_ */
