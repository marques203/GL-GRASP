/*
 * GRASPv2.h
 *
 *  Created on: 23 mag 2017
 *      Author: antonio
 */

#include "GRASP.h"

#ifndef GRASPV2_H_
#define GRASPV2_H_

class GRASP_v2 : public GRASP {
public:
	GRASP_v2(HDAG &I, char *argv[]);
	virtual ~GRASP_v2();

	virtual void algorithm(HDAG &I, int argc, char *argv[]);
	virtual void run(HDAG &I, HDAG &S);
	virtual void construction(HDAG &I, HDAG &S);
	virtual void allocate_solution(HDAG &I, HDAG &S);
	virtual void add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);
	virtual unsigned check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	virtual unsigned find_feasible_position_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);


};

#endif /* GRASPV2_H_ */
