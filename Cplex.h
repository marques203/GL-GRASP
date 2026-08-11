/*
 * Cplex.h
 *
 *  Created on: 26 mag 2017
 *      Author: antonio
 */

#ifndef CPLEX_H_
#define CPLEX_H_

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "HDAG.h"
#include <string.h>
#include "SolutionIMLCM.h"
#include "chrono"

using namespace std;

class Cplex {
public:
	Cplex(char *argv[]);
	virtual ~Cplex();

	void cplex(HDAG &G, int argc, char *argv[]);

	void solve_model(HDAG &G);

	void build_edge_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_node_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_position_structure(HDAG &G, std::vector<unsigned> &lev_gap,
			IloIntArray &P);

	// c^t_u1v1u2v2 <= x^{t+1}_v1v2 - x^t_u1u2 <= - c^t_u1v1u2v2
	void add_first_constraints(unsigned u1, unsigned v1, unsigned u2,
			unsigned v2, IloBoolVarArray &C, IloBoolVarArray &X, IloModel &model,
			unsigned i, SolutionIMLCM &S);

	// 1 - c^t_u1v1u2v2 <= x^{t+1}_v2v1 + x^t_u1u2 <= 1 + c^t_u1v1u2v2
	void add_second_constraints(unsigned u1, unsigned v1, unsigned u2,
			unsigned v2, IloBoolVarArray &C, IloBoolVarArray &X, IloModel &model,
			unsigned i, SolutionIMLCM &S);

	// 0 <= X^t_{u1u2} + X^t{u2u3} - X^t_{u1u3} <= 1
	void add_third_constraints(unsigned l, unsigned i, unsigned j, unsigned k,
			IloBoolVarArray &X, IloModel &model, HDAG &G,
			std::vector<unsigned> &lev_gap, SolutionIMLCM &S);

	void add_incremental_constraints(SolutionIMLCM &S, IloBoolVarArray &X,
			IloModel &model, HDAG &G, std::vector<unsigned> &lev_gap);

	void add_position_constraints(SolutionIMLCM &S, HDAG &G,
			std::vector<unsigned> &lev_gap, IloModel &model, IloIntArray &P,
			IloIntExprArray &P_, IloBoolVarArray &X, IloEnv &env);

	bool edge_pair_already_selected(unsigned i, unsigned j, unsigned k,
			unsigned l);

	void add_edge_pair(unsigned l, unsigned i, unsigned j, unsigned k,
			IloBoolVarArray &X, IloModel &model, HDAG &G,
			std::vector<unsigned> &lev_gap, SolutionIMLCM &S);

	bool node_pair_already_selected(unsigned i, unsigned j);

	void add_node_pair(unsigned i, unsigned j);

	unsigned find_node_pair(unsigned i, unsigned j);

	bool optimal_solution_found() { return optimal; };

	double get_runtime() { return runtime; };
	void set_runtime(double x) { runtime = x; };


private:
	unsigned k;
	int tl;
	double runtime;
	bool optimal;

};

#endif /* CPLEX_H_ */
