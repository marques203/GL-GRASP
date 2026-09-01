/*
 * Localsolver.h
 *
 *  Created on: 28 set 2017
 *      Author: antonio
 */

#ifndef LOCALSOLVER_H_
#define LOCALSOLVER_H_

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "HDAG.h"
#include <string.h>
#include "SolutionIMLCM.h"
#include "chrono"
#include "localsolver.h"

using namespace localsolver;
using namespace std;

class Localsolver {
public:
	Localsolver(char *argv[]);
	virtual ~Localsolver();

	void localsolver(HDAG &G, int argc, char *argv[]);

	void solve_model(HDAG &G, char *outfile);

	void build_edge_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_node_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_position_structure(HDAG &G, std::vector<unsigned> &lev_gap,
			vector<lsint> &P);

	// c^t_u1v1u2v2 <= x^{t+1}_v1v2 - x^t_u1u2 <= - c^t_u1v1u2v2
	void add_first_constraints(unsigned u1, unsigned v1, unsigned u2,
			unsigned v2, vector<localsolver::LSExpression> &C, vector<localsolver::LSExpression> &X,
			LSModel &model, unsigned i, SolutionIMLCM &S);

	// 1 - c^t_u1v1u2v2 <= x^{t+1}_v2v1 + x^t_u1u2 <= 1 + c^t_u1v1u2v2
	void add_second_constraints(unsigned u1, unsigned v1, unsigned u2,
			unsigned v2, vector<localsolver::LSExpression> &C, vector<localsolver::LSExpression> &X,
			LSModel &model, unsigned i, SolutionIMLCM &S);

	// 0 <= X^t_{u1u2} + X^t{u2u3} - X^t_{u1u3} <= 1
	void add_third_constraints(unsigned l, unsigned i, unsigned j, unsigned k,
			vector<localsolver::LSExpression> &X, LSModel &model, HDAG &G,
			std::vector<unsigned> &lev_gap, SolutionIMLCM &S);

	void add_incremental_constraints(SolutionIMLCM &S, vector<localsolver::LSExpression> &X,
			LSModel &model, HDAG &G, std::vector<unsigned> &lev_gap);

	void add_position_constraints(SolutionIMLCM &S, HDAG &G,
			std::vector<unsigned> &lev_gap, LSModel &model, vector<lsint> &P,
			vector<localsolver::LSExpression> &P_, vector<localsolver::LSExpression> &X);

	double get_runtime() { return runtime; };
	void set_runtime(double x) { runtime = x; };


private:
	unsigned k;
	int tl;
	double runtime;


};

#endif /* LOCALSOLVER_H_ */
