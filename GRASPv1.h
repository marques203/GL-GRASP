/*
 * GRASPv1.h
 *
 *  Created on: 24 mag 2017
 *      Author: antonio
 */

#include "GRASP.h"

#ifndef GRASPV1_H_
#define GRASPV1_H_

class GRASP_v1 : public GRASP {
public:
	GRASP_v1(HDAG &I, char *argv[]);
	virtual ~GRASP_v1();

	virtual void algorithm(HDAG &I, int argc, char *argv[]);
	virtual void run(HDAG &I, HDAG &S);
	virtual void construction(HDAG &I, HDAG &S);
	virtual void allocate_solution(HDAG &I, HDAG &S);
	virtual void add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);
	virtual unsigned check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	virtual unsigned find_feasible_position_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	virtual unsigned find_feasible_position_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	virtual unsigned find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);

	void add_first_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);
	bool range_feasibility_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	bool range_feasibility_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);
	unsigned find_free_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);

private:

	std::vector<unsigned> upmost_original;
	std::vector<unsigned> original_to_assign;


};

#endif /* GRASPV1_H_ */
