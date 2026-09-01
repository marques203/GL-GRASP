/*
 * TABU.h
 *
 *  Created on: 19 lug 2017
 *      Author: antonio
 */

#include "GRASP.h"
#include <tuple>

#ifndef TABU_H_
#define TABU_H_

using namespace std;

class TABU : public GRASP {
public:
	TABU(HDAG &I, char *argv[]);
	virtual ~TABU();

	virtual void algorithm(HDAG &I, int argc, char *argv[]);
	virtual void run(HDAG &I, HDAG &S) {};
	void run_TABU(HDAG &I, HDAG &S, unsigned iter);
	virtual void construction(HDAG &I, HDAG &S);
	virtual void allocate_solution(HDAG &I, HDAG &S);
	virtual void add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng) {};
	virtual void add_node_to_solution_for_TABU(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned new_pos);
	virtual unsigned check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned pos);
	virtual unsigned find_feasible_position_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };

	void build_UNASSIGNED(HDAG &I, HDAG &S);
	double compute_Bward_potential_cost(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned new_p);
	double compute_Fward_potential_cost(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned new_p);
	unsigned find_new_position(unsigned l, unsigned i);
	void build_CL_TABU(HDAG &I, HDAG &S);
	pair<unsigned, unsigned> extract_from_CL_TABU();
	void apply_updates_after_insertion(HDAG &I, HDAG &S, unsigned l, unsigned i);
	void update_CL_TABU(HDAG &I, HDAG &S);

	void update_frequencies_matrix(HDAG &I, HDAG &S);

	void free_UNASSIGNED() {
		for (unsigned l = 0; l < UNASSIGNED.size(); ++l)
			for (unsigned i = 0; i < UNASSIGNED[l].size(); ++i)
				while (UNASSIGNED[l][i].size() > 0)
					UNASSIGNED[l][i].pop_back();
		return;
	}

	double getTau() { return tau; }
	void updateTau() {
		tau = global_min_g + getalpha()*(global_max_g - global_min_g);
	} ;


	void build_RCL_for_TABU(HDAG &I, HDAG &S);
	void re_build_RCL_for_TABU(HDAG &I, HDAG &S, double old_tau, unsigned l, unsigned i);

private:

	double beta;
	double tau;
	vector<vector<vector<tuple<unsigned, unsigned, double, bool>>>> UNASSIGNED;
	double global_min_g;
	double global_max_g;
	vector<vector<double>> min_g;
	tuple<unsigned, unsigned, unsigned> position_min_g;
	tuple<unsigned, unsigned, unsigned> position_max_g;
	vector<vector<vector<unsigned>>> FREQUENCIES;
	double min_value_in_RCL;
	pair<unsigned, unsigned> minimum_in_RCL;

};

#endif /* TABU_H_ */
