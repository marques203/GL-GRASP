/*
 * Cplex.cpp
 *
 *  Created on: 26 mag 2017
 *      Author: antonio
 */

#include "Cplex.h"
#include <ilconcert/iloalg.h>
#include <ilconcert/iloenv.h>
#include <ilconcert/iloexpression.h>
#include <ilconcert/iloextractable.h>
#include <ilconcert/ilolinear.h>
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <ilcplex/ilocplexi.h>
#include <vector>

Cplex::Cplex(char *argv[]): k(0), tl(100), runtime(0), optimal(false) {
	// TODO Auto-generated constructor stub

	k = atoi(argv[3]);
	tl = atoi(argv[6]);

}

Cplex::~Cplex() {
	// TODO Auto-generated destructor stub
}

void Cplex::cplex(HDAG &G, int argc, char *argv[]) {

	SolutionIMLCM S;
	double initial_cost = S.getCost(G);

	char *outfile = argv[5];
	FILE *out_f;
	out_f = fopen(outfile, "a");
	auto start_chrono = chrono::system_clock::now();
	solve_model(G);
	auto end_chrono = chrono::system_clock::now();
	auto elapsed = chrono::duration_cast<chrono::milliseconds>(end_chrono - start_chrono);
	set_runtime((double)((double)elapsed.count() / (double)1000));
	if (optimal_solution_found()) {
		if (get_runtime() > tl + 100) {
			fprintf(out_f, "%0.3f %0.3f %d\n", initial_cost, (double)(tl), 0);
		}
		else {
			fprintf(out_f, "%0.3f %0.3f %d\n", G.get_cost(), get_runtime(), 1);
		}
	}
	else {
		if (get_runtime() > tl + 100) {
			fprintf(out_f, "%0.3f %0.3f %d\n", initial_cost, (double)(tl), 0);
		}
		else {
			fprintf(out_f, "%0.3f %0.3f %d\n", G.get_cost(), get_runtime(), 0);
		}
	}
	fclose(out_f);

//	G.printHDAG(G, "HDAG-cplex.tex");

}

void Cplex::solve_model(HDAG &G) {

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

	IloEnv env;
	IloModel model(env);

	IloBoolVarArray X(env, N.size());
	IloBoolVarArray C(env, M.size());
	IloIntArray P(env, n_original_nodes); // Initial position of the original nodes
	IloIntExprArray P_(env, n_original_nodes);	// Final position of the original nodes

	build_position_structure(G, lev_gap, P);

	for (unsigned i = 0; i < N.size(); ++i) {
		X[i] = IloBoolVar(env);
	}

	for (unsigned i = 0; i < M.size(); ++i) {
		C[i] = IloBoolVar(env);
	}

	for (unsigned i = 0; i < n_original_nodes; ++i) {
		P_[i] = IloIntExpr(env);
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

	add_position_constraints(S, G, lev_gap, model, P, P_, X, env);

	IloIntExpr sum_over_all(env, 0);

	for (unsigned i = 0; i < M.size(); ++i) {
		sum_over_all += C[i];
	}

	model.add(IloMinimize(env, sum_over_all));

	IloCplex cplex(model);

	cplex.setParam(IloCplex::TiLim, tl);

	cplex.setOut(env.getNullStream());

	cplex.solve();

	IloCplex::Status st = cplex.getCplexStatus();

//	cout << " status = " << st << endl;
//	cout << " min (cplex) = " << cplex.getObjValue() << endl;

//	for (unsigned elem = 0; elem < n_original_nodes; ++elem)
//		cout << "P=" << P[elem] << "(" << cplex.getValue(P_[elem]) << ")" << endl;

	if (st != IloCplex::AbortTimeLim) {

		optimal = true;

	}
	else {

		set_runtime(tl);

	}


	G.set_cost(cplex.getObjValue());

	unsigned change = N.size();

	while (change) {

		change = 0;

		for (unsigned i = 0; i < N.size(); ++i) {

			pair<unsigned, unsigned> pos_i =
					S.rescale_to_original_and_get_position(G, lev_gap,
							N[i].first);
			pair<unsigned, unsigned> pos_j =
					S.rescale_to_original_and_get_position(G, lev_gap,
							N[i].second);

			if (cplex.getValue(X[i]) == 0) {

				if (pos_i.second < pos_j.second) {

					G.swapPositions(pos_i.second, pos_j.second, pos_i.first);

					change++;

				}

			}
			if (cplex.getValue(X[i]) == 1) {

				if (pos_i.second > pos_j.second) {

					G.swapPositions(pos_i.second, pos_j.second, pos_i.first);

					change++;

				}

			}

		}

	}

	double cost = 0;
	for (unsigned i = 0; i < M.size(); ++i) {
		if (round(cplex.getValue(C[i])) == 1)
			cost++;
	}

/*	cout << "cost (our) = " << G.get_cost() << endl;
	cout << "cost (manual) = " << cost << endl;*/


/*	for (unsigned i = 0; i < n_original_nodes; ++i) {
		cout << " old pos = " << P[i] << " - new pos = " << cplex.getValue(P_[i])<< endl;

	}*/

/*	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size() - 1; ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				for (unsigned i_ = i + 1; i_ < LEVELS[l].size(); ++i_) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// arc (i, j)
						unsigned id_i = IDs[l][i];
						unsigned id_j = LEVELS[l][i][j];

						// arc (i_, j_)
						unsigned id_i_ = IDs[l][i_];
						unsigned id_j_ = LEVELS[l][i_][j_];

						if (G.areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

							for (unsigned x = 0; x < M.size(); ++x) {

								pair<unsigned, unsigned > i_p;
								pair<unsigned, unsigned > j_p;
								pair<unsigned, unsigned > k_p;
								pair<unsigned, unsigned > l_p;
								i_p = S.rescale_to_original_and_get_position(G, lev_gap, M[x].first.first);
								j_p = S.rescale_to_original_and_get_position(G, lev_gap, M[x].first.second);
								k_p = S.rescale_to_original_and_get_position(G, lev_gap, M[x].second.first);
								l_p = S.rescale_to_original_and_get_position(G, lev_gap, M[x].second.second);

								unsigned id_i_p = IDs[i_p.first][i_p.second];
								unsigned id_j_p = IDs[j_p.first][j_p.second];
								unsigned id_k_p = IDs[k_p.first][k_p.second];
								unsigned id_l_p = IDs[l_p.first][l_p.second];

								if (id_i == id_i_p && id_j == id_j_p &&
									id_i_ == id_k_p && id_j_ == id_l_p &&
									l == i_p.first) {

									if (cplex.getValue(C[x]) == 0) {

									   cout << " in level " << l << endl;
									   cout << " pair (" << id_i << ", " << id_j << ") - (" <<
											id_i_ << ", " << id_j_ << ") not present!" << endl;

									}

									break;

								}

							}

						}

					}

				}

			}

		}

	}

	int count = 0;

	for (unsigned x = 0; x < M.size(); ++x) {

		pair<unsigned, unsigned > i;
		pair<unsigned, unsigned > j;
		pair<unsigned, unsigned > k;
		pair<unsigned, unsigned > l;
		i = S.rescale_to_original_and_get_position(G, lev_gap, M[x].first.first);
		j = S.rescale_to_original_and_get_position(G, lev_gap, M[x].first.second);
		k = S.rescale_to_original_and_get_position(G, lev_gap, M[x].second.first);
		l = S.rescale_to_original_and_get_position(G, lev_gap, M[x].second.second);

		unsigned id_i = IDs[i.first][i.second];
		unsigned id_j = IDs[j.first][j.second];
		unsigned id_k = IDs[k.first][k.second];
		unsigned id_l = IDs[l.first][l.second];


		if (G.areCrossingEdge(id_i, id_j, id_k, id_l, i.first)) {
			if (cplex.getValue(C[x])) {
		//		cout << " CONGRUENT!" << endl;
			}
			else {
				cout << " ARE CROSSING!!! <----------------------------------------------" << endl;
				cout << " PAIR of edges " << x << ": ("
						<< id_i << ", " << id_j << ") - ("
						<< id_k << ", " << id_l << ") o ("
						<< M[x].first.first << ", " << M[x].first.second << ") - ("
						<< M[x].second.first << ", " << M[x].second.second << ") ";
				cout << " - variable: " << C[x].getName() << " = " << cplex.getValue(C[x]);
				cout << " CROSS: y - ";
				cout << " INCONGRUENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
			}
		}
		else {
			if (!cplex.getValue(C[x])) {
		//		cout << " CONGRUENT!" << endl;
			}
			else {
				cout << " NOT ARE CROSSING!!!" << endl;
				cout << " PAIR of edges " << x << ": ("
						<< id_i << ", " << id_j << ") - ("
						<< id_k << ", " << id_l << ") o ("
						<< M[x].first.first << ", " << M[x].first.second << ") - ("
						<< M[x].second.first << ", " << M[x].second.second << ") ";
				cout << " - variable: " << C[x].getName() << " = " << cplex.getValue(C[x]);
				cout << " CROSS: n - ";
				cout << " INCONGRUENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! " << endl;
				count ++;
			}
		}

	}

	cout << " total incongruents " << count << endl;

	*/

	env.end();

}

void Cplex::build_position_structure(HDAG &G, vector<unsigned> &lev_gap,
		IloIntArray &P) {

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

void Cplex::add_first_constraints(unsigned u1, unsigned v1, unsigned u2,
		unsigned v2, IloBoolVarArray &C, IloBoolVarArray &X, IloModel &model,
		unsigned i, SolutionIMLCM &S) {

	unsigned p_1 = S.find_node_pair(v1, v2);

	unsigned p_2 = S.find_node_pair(u1, u2);

	model.add(X[p_1] - X[p_2] <= C[i]);
	model.add(X[p_1] - X[p_2] >= -C[i]);

	string name = "c" + to_string(u1) + to_string(v1) + to_string(u2)
			+ to_string(v2);


	C[i].setName(name.c_str());

	name = "";
	name = "x" + to_string(p_1);
	X[p_1].setName(name.c_str());

	name = "";
	name = "x" + to_string(p_2);
	X[p_2].setName(name.c_str());

}

void Cplex::add_second_constraints(unsigned u1, unsigned v1, unsigned u2,
		unsigned v2, IloBoolVarArray &C, IloBoolVarArray &X, IloModel &model,
		unsigned i, SolutionIMLCM &S) {

	unsigned p_1 = S.find_node_pair(v2, v1);

	unsigned p_2 = S.find_node_pair(u1, u2);

	model.add(X[p_1] + X[p_2] <= 1 + C[i]);
	model.add(X[p_1] + X[p_2] >= 1 - C[i]);

//	model.add(1 - C[i] <= (X[p_1] + X[p_2]) <= 1 + C[i]);

	string name = "c" + to_string(u1) + to_string(v1) + to_string(u2)
			+ to_string(v2);

	C[i].setName(name.c_str());

	name = "";
	name = "x" + to_string(p_1);
	X[p_1].setName(name.c_str());

	name = "";
	name = "x" + to_string(p_2);
	X[p_2].setName(name.c_str());

}

void Cplex::add_third_constraints(unsigned l, unsigned i, unsigned j, unsigned k,
		IloBoolVarArray &X, IloModel &model, HDAG &G,
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

//		unsigned p_1 = S.find_node_pair(u1, u2);
//		unsigned p_2 = S.find_node_pair(u2, u3);
//		unsigned p_3 = S.find_node_pair(u1, u3);

		// 0 <= X^t_{u1u2} + X^t{u2u3} - X^t_{u1u3} <= 1
		model.add(X[p_1] + X[p_2] - X[p_3] <= 1);
		model.add(X[p_1] + X[p_2] - X[p_3] >= 0);

//		model.add(0 <= (X[p_1] + X[p_2]) - X[p_3] <= 1);

		string name;
		name = "";
		name = "x" + to_string(u1) + to_string(u2);
		X[p_1].setName(name.c_str());

		name = "";
		name = "x" + to_string(u2) + to_string(u3);
		X[p_2].setName(name.c_str());

		name = "";
		name = "x" + to_string(u1) + to_string(u3);
		X[p_3].setName(name.c_str());

	}

	return;

}

void Cplex::add_incremental_constraints(SolutionIMLCM &S,
		IloBoolVarArray &X, IloModel &model, HDAG &G,
		std::vector<unsigned> &lev_gap) {

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
					model.add(X[p] == 1);
				else
					//if (u1.second > u2.second)
					model.add(X[p] == 0);
			}
		}

	}

}

void Cplex::add_position_constraints(SolutionIMLCM &S, HDAG &G,
		std::vector<unsigned> &lev_gap, IloModel &model, IloIntArray &P,
		IloIntExprArray &P_, IloBoolVarArray &X, IloEnv &env) {

	auto &LEVELS = G.getLEVELS();
	auto &N = S.getN();
	IloExpr prova(env);

	unsigned elem = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {

		unsigned nl = 0;
		for (unsigned i = 0; i < LEVELS[l].size(); ++i)
				nl++;

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			if (G.isOriginalNode(i, l)) {

				P_[elem] = P_[elem] + ((int)nl - 1);

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

										P_[elem] = P_[elem] - X[x];

										break;

									}

									if (p1.second == i_ && p2.second == i) {

										P_[elem] = P_[elem] - (1 - X[x]);

										break;

									}

								}

							}

					}

				}

				int n_ = LEVELS[l].size() - 1;
				IloInt lb = IloMax(0, P[elem] - k);
				IloInt ub = IloMin(P[elem] + k, n_);

//				model.add(lb <= P_[elem]
//							<= ub);

				model.add(P_[elem] >= lb);
				model.add(P_[elem] <= ub);

				elem++;

			}

		}

	}

}
