/*
 * Localsolver.cpp
 *
 *  Created on: 28 set 2017
 *      Author: antonio
 */

#include "Localsolver.h"

using namespace localsolver;
using namespace std;

Localsolver::Localsolver(char *argv[]): k(0), tl(100), runtime(0) {
	// TODO Auto-generated constructor stub

	k = atoi(argv[3]);
	tl = atoi(argv[6]);

}

Localsolver::~Localsolver() {
	// TODO Auto-generated destructor stub
}

void Localsolver::localsolver(HDAG &G, int argc, char *argv[]) {

	solve_model(G, argv[5]);

}

void Localsolver::solve_model(HDAG &G, char *outfile) {

	SolutionIMLCM s;
	double cost = s.getCost(G);

	auto start_chrono = chrono::system_clock::now();

	FILE *out_f;
	out_f = fopen(outfile, "a");

	SolutionIMLCM S;
	G.set_cost(S.getCost(G));

	auto &LEVELS = G.getLEVELS();
	auto &IDs = G.getIDs();
	const auto &N = S.getN();
	const auto &M = S.getM();
	unsigned n_original_nodes = G.get_n_original_nodes(G);

	vector<unsigned> lev_gap = S.compute_level_gap_vector(G);

	S.build_edge_pair_structure(G, lev_gap);

	S.build_node_pair_structure(G, lev_gap);

    LocalSolver localsolver;
    LSModel model = localsolver.getModel();

    vector<localsolver::LSExpression> X;  X.resize(N.size());
    vector<localsolver::LSExpression> C;  C.resize(M.size());
    vector<lsint> P;					  P.resize(n_original_nodes);
    vector<localsolver::LSExpression> P_; P_.resize(n_original_nodes);

    build_position_structure(G, lev_gap, P);


    for (unsigned i = 0; i < N.size(); ++i) {
    	X[i] = model.boolVar();
    }


    for (unsigned i = 0; i < M.size(); ++i) {
    	C[i] = model.boolVar();
    }

	for (unsigned i = 0; i < M.size(); i++) {

		unsigned u1 = M[i].first.first;
		unsigned v1 = M[i].first.second;
		unsigned u2 = M[i].second.first;
		unsigned v2 = M[i].second.second;

		if (v1 < v2) {

			add_first_constraints(u1, v1, u2, v2, C, X, model, i, S);

		}

		if (v1 > v2) {

			add_second_constraints(u1, v1, u2, v2, C, X, model, i, S);

		}

	}

	for (unsigned l = 0; l < G.getLevNumber(); ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l].size(); ++j) {

				for (unsigned k = 0; k < LEVELS[l].size(); ++k) {

					add_third_constraints(l, i , j, k, X, model, G, lev_gap, S);

				}

			}

		}
	}

	add_incremental_constraints(S, X, model, G, lev_gap);

	add_position_constraints(S, G, lev_gap, model, P, P_, X);

	LSExpression sum_over_all = model.sum();
	sum_over_all += 0;
	for (unsigned i = 0; i < M.size(); ++i)
		sum_over_all += C[i];

	model.minimize(sum_over_all);

	model.close();

    LSPhase phase = localsolver.createPhase();
    phase.setTimeLimit(tl);

	localsolver.solve();

	if ( 2 == localsolver.getSolution().getStatus() )
		cost = (double)(sum_over_all.getValue());

//	for (unsigned elem = 0; elem < n_original_nodes; ++elem)
//		cout << "P=" << P[elem] << "(" << P_[elem].getValue() << ")" << endl;

//	for (unsigned i = 0; i < N.size(); ++i)
//		cout << "X= " << X[i].getValue() << endl;

	auto end_chrono = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>(end_chrono - start_chrono);
	set_runtime((double)((double)elapsed.count() / (double)1000));
	fprintf(out_f, "%0.3f %0.3f %d\n", cost, get_runtime(), 0);
	fclose(out_f);

	return;

}

void Localsolver::build_position_structure(HDAG &G, vector<unsigned> &lev_gap,
		vector<lsint> &P) {

	auto &LEVELS = G.getLEVELS();
	unsigned elem = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		unsigned n_l = 0;
		for (unsigned i = 0; i < LEVELS[l].size(); ++i)
			if (G.isOriginalNode(i, l))
				n_l++;

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			if (G.isOriginalNode(i, l)) {
				P[elem] = n_l;
				for (unsigned i_ = i; i_ < LEVELS[l].size(); ++i_) {
					if (G.isOriginalNode(i_, l)) {
						P[elem] = P[elem] - 1;
					}
				}
				elem++;
			}
		}
	}
}


void Localsolver::add_first_constraints(unsigned u1, unsigned v1, unsigned u2,
		unsigned v2, vector<localsolver::LSExpression> &C, vector<localsolver::LSExpression> &X,
		LSModel &model, unsigned i, SolutionIMLCM &S) {

	unsigned p_1 = S.find_node_pair(v1, v2);

	unsigned p_2 = S.find_node_pair(u1, u2);

	model.addConstraint(X[p_1] - X[p_2] <= 0 + C[i]);

	model.addConstraint(X[p_1] - X[p_2] >= 0 - C[i]);

}

void Localsolver::add_second_constraints(unsigned u1, unsigned v1, unsigned u2,
		unsigned v2, vector<localsolver::LSExpression> &C, vector<localsolver::LSExpression> &X,
		LSModel &model, unsigned i, SolutionIMLCM &S) {

	unsigned p_1 = S.find_node_pair(v2, v1);

	unsigned p_2 = S.find_node_pair(u1, u2);

	model.addConstraint(X[p_1] + X[p_2] <= 1 + C[i]);
	model.addConstraint(X[p_1] + X[p_2] >= 1 - C[i]);

}


void Localsolver::add_third_constraints(unsigned l, unsigned i, unsigned j, unsigned k,
			vector<localsolver::LSExpression> &X, LSModel &model, HDAG &G,
			std::vector<unsigned> &lev_gap, SolutionIMLCM &S) {

	auto &N = S.getN();
	unsigned u1 = G.getIDs()[l][i] + lev_gap[l];
	unsigned u2 = G.getIDs()[l][j] + lev_gap[l];
	unsigned u3 = G.getIDs()[l][k] + lev_gap[l];

	unsigned p_1 = 0;
	unsigned p_2 = 0;
	unsigned p_3 = 0;

	bool u1_u2_found = false;
	bool u2_u3_found = false;
	bool u1_u3_found = false;
	bool all_found = false;

	for (unsigned x = 0; x < N.size(); ++x) {
		if (N[x].first == u1 && N[x].second == u2) {
			u1_u2_found = true;
			p_1 = x;
		}
		if (N[x].first == u2 && N[x].second == u3) {
			u2_u3_found = true;
			p_2 = x;
		}
		if (N[x].first == u1 && N[x].second == u3) {
			u1_u3_found = true;
			p_3 = x;
		}
		if (u1_u2_found && u2_u3_found && u1_u3_found) {
			all_found = true;
			break;
		}
	}

	if (!all_found) return;

	if (u1 < u2 && u2 < u3) {

		// 0 <= X^t_{u1u2} + X^t{u2u3} - X^t_{u1u3} <= 1
		model.addConstraint(X[p_1] + X[p_2] - X[p_3] <= 1);
		model.addConstraint(X[p_1] + X[p_2] - X[p_3] >= 0);

	}

	return;

}

void Localsolver::add_incremental_constraints(SolutionIMLCM &S, vector<localsolver::LSExpression> &X,
		LSModel &model, HDAG &G, std::vector<unsigned> &lev_gap) {

	const auto &N = S.getN();
	const auto &IDs = G.getIDs();

	for (unsigned i = 0; i < N.size(); ++i) {

		pair<unsigned, unsigned> u1 = S.rescale_to_original_and_get_position(G,
				lev_gap, N[i].first);
		pair<unsigned, unsigned> u2 = S.rescale_to_original_and_get_position(G,
				lev_gap, N[i].second);

		if (G.isOriginalNode(IDs[u1.first][u1.second], u1.first)) {
			if (G.isOriginalNode(IDs[u2.first][u2.second], u2.first)) {
				unsigned p = S.find_node_pair(N[i].first, N[i].second);
				if (u1.second < u2.second)
					model.addConstraint(X[p] == 1);
				else
					//if (u1.second > u2.second)
					model.addConstraint(X[p] == 0);
			}
		}

	}

}


void Localsolver::add_position_constraints(SolutionIMLCM &S, HDAG &G,
		std::vector<unsigned> &lev_gap, LSModel &model, vector<lsint> &P,
		vector<localsolver::LSExpression> &P_, vector<localsolver::LSExpression> &X) {

	auto &LEVELS = G.getLEVELS();
	auto &N = S.getN();

	unsigned elem = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {

		unsigned nl = 0;
		for (unsigned i = 0; i < LEVELS[l].size(); ++i)
				nl++;

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			if (G.isOriginalNode(i, l)) {

				lsint lb;
				if (0 > P[elem] - k) {
					lb = 0;
				}
				else {
					lb = P[elem] - k;
				}

				lsint ub;
				if (P[elem] + k < nl) {
					ub = P[elem] + k;
				}
				else {
					ub = nl;
				}

				LSExpression quantity = model.sum();
				quantity += 0;
				quantity += (int(nl) - 1);

//				P_[elem] += 0;
//				P_[elem] += (int(nl) - 1);

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); ++i_) {

					if (i != i_) {

							for (unsigned x = 0; x < N.size(); ++x) {

								pair<unsigned, unsigned> p1 =
										S.rescale_to_original_and_get_position(
												G, lev_gap, N[x].first);
								pair<unsigned, unsigned> p2 =
										S.rescale_to_original_and_get_position(
												G, lev_gap, N[x].second);

								if ((p1.first == p2.first) && (p2.first == l)) {

									if (p1.second == i && p2.second == i_) {

										//P_[elem].operator-(X[x]);

//										P_[elem] += (0 - X[x]);
//										P_[elem] += not(X[x]);

										quantity = quantity - X[x];

										break;

									}

									if (p1.second == i_ && p2.second == i) {

										//P_[elem].operator-(1 - X[x]);

//										P_[elem] += (0 - (1 - X[x]));

//										P_[elem] += (1 + not(X[x]));

										quantity = quantity - (1 - X[x]);

										break;

									}

								}

							}

					}

				}

				P_[elem] = model.intVar(lb, ub);

				model.addConstraint(P_[elem] == quantity);
				model.addConstraint(P_[elem] >= lb);
				model.addConstraint(P_[elem] <= ub);

				elem++;

			}

		}

	}

}


