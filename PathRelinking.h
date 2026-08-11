/*
 * PathRelinking.h
 *
 *  Created on: 12 giu 2017
 *      Author: antonio
 */

#include "EliteSet.h"
#include "SolutionIMLCM.h"
#include <string.h>
#include "chrono"

#ifndef PATHRELINKING_H_
#define PATHRELINKING_H_

using namespace std;

class Path_Relinking {
public:
	Path_Relinking(char *c);
	virtual ~Path_Relinking();

	void PR(Elite_Set &ES, HDAG &best_solution);

	void Relink(Elite_Set &ES, HDAG &best_solution, HDAG &I);

	void forward(HDAG &S1, HDAG &S2, HDAG &best, HDAG &I);
	void backward(HDAG &S1, HDAG &S2, HDAG &best);
	void mixed(HDAG &S1, HDAG &S2, HDAG &best);

	void forward_step(HDAG &source, HDAG &target, HDAG &best, HDAG &I);

	void subtract_cost_in_level(HDAG &S, unsigned l, double &cost);
	void add_cost_in_level(HDAG &S, unsigned l, double &cost);

	char *get_PR_type() { return PR_type; };

private:

	char *PR_type;

};

#endif /* PATHRELINKING_H_ */
