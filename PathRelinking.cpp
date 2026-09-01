/*
 * PathRelinking.cpp
 *
 *  Created on: 12 giu 2017
 *      Author: antonio
 */

#include "PathRelinking.h"
#include "limits.h"

Path_Relinking::Path_Relinking(char *c) : PR_type() {
	// TODO Auto-generated constructor stub

	PR_type=c;

}

void Path_Relinking::Relink(Elite_Set &ES, HDAG &best_solution, HDAG &I) {

	auto &EliteSet = ES.getElite();

	HDAG *solution = new HDAG;
	solution->copy_HDAG(*EliteSet[0]);
	solution->set_cost(INT_MAX);
//	double old_best_cost_dummy = solution->get_cost();
//	double old_overall_best_dummy= best_solution.get_cost();

	if (strcmp(PR_type, "f") == 0 || strcmp(PR_type, "f") == 0) {
		for (unsigned i = 0; i < EliteSet.size() - 1; ++i) {
			for (unsigned i_ = i + 1; i_ < EliteSet.size(); ++i_) {

/*				if (EliteSet[i]->get_cost() > EliteSet[i_]->get_cost()) {
					cout << " Relink " << "(" << EliteSet[i]->get_cost() <<
							") --f-> (" << EliteSet[i_]->get_cost() << ")";

				}
				else {
					cout << " Relink " << "(" << EliteSet[i_]->get_cost() <<
							") --f-> (" << EliteSet[i]->get_cost() << ")";

				}*/

				forward(*EliteSet[i], *EliteSet[i_], *solution, I);

				// controllo con best_solution


/*				if (old_best_cost_dummy > solution->get_cost()) {

					if (old_overall_best_dummy > solution->get_cost()) {

						cout << " ---> <" << solution->get_cost() << "> "
								<< " ---> *new best* " << endl;
						old_overall_best_dummy = solution->get_cost();

					}
					else {
						cout << " ---> <" << solution->get_cost() << ">" << endl;
					}
					old_best_cost_dummy = solution->get_cost();
				}
				else {
					cout << endl;
				}*/

				if (solution->get_cost() < best_solution.get_cost()) {
					for (unsigned l = 0; l < solution->getLEVELS().size(); ++l) {
						best_solution.copy_level(*solution, l);
					}

					best_solution.set_cost(solution->get_cost());
					FILE *FF = fopen("Relinking", "a");
					fprintf(FF, "%f\n", best_solution.get_cost());
					fclose(FF);

				}

			}
		}
	}

/*	if (best_solution.get_cost() > solution->get_cost()) {

		for (unsigned l = 0; l < solution->getLEVELS().size(); ++l) {
			best_solution.copy_level(*solution, l);
		}
		best_solution.set_cost(solution->get_cost());

		solution->printHDAG(*solution, "HDAG-tabu+PR.tex");

	}*/

}

void Path_Relinking::forward(HDAG &S1, HDAG &S2, HDAG &best, HDAG &I) {

	if (S1.get_cost() > S2.get_cost())
		forward_step(S1, S2, best, I);
	else
		forward_step(S2, S1, best, I);

	return;

}

void Path_Relinking::subtract_cost_in_level(HDAG &S, unsigned l, double &cost) {


	const auto &Lev_S = S.getLEVELS();

	for (unsigned i = 0; i < Lev_S[l].size() - 1; ++i) {

		for (unsigned j = 0; j < Lev_S[l][i].size(); ++j) {

			for (unsigned i_ = i + 1; i_ < Lev_S[l].size(); ++i_) {

				for (unsigned j_ = 0; j_ < Lev_S[l][i_].size(); ++j_) {

					// Arc (i, j)
					unsigned id_i = S.getIDs()[l][i];
					unsigned id_j = Lev_S[l][i][j];


					// Arc (i_, j_)
					unsigned id_i_ = S.getIDs()[l][i_];
					unsigned id_j_ = Lev_S[l][i_][j_];


					if (S.areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {

						cost--;

					}
				}
			}
		}
	}

	return;

}

void Path_Relinking::add_cost_in_level(HDAG &S, unsigned l, double &cost) {

	const auto &Lev_S = S.getLEVELS();

	for (unsigned i = 0; i < Lev_S[l].size() - 1; ++i) {

		for (unsigned j = 0; j < Lev_S[l][i].size(); ++j) {

			for (unsigned i_ = i + 1; i_ < Lev_S[l].size(); ++i_) {

				for (unsigned j_ = 0; j_ < Lev_S[l][i_].size(); ++j_) {

					// Arc (i, j)
					unsigned id_i = S.getIDs()[l][i];
					unsigned id_j = Lev_S[l][i][j];

					// Arc (i_, j_)
					unsigned id_i_ = S.getIDs()[l][i_];
					unsigned id_j_ = Lev_S[l][i_][j_];

					if (S.areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
						cost++;
					}
				}
			}
		}
	}

	return;

}

void Path_Relinking::forward_step(HDAG &source, HDAG &target, HDAG &best, HDAG &I) {

/*	double actual_cost;
	HDAG *S = new HDAG;
	auto &LEVELS = target.getLEVELS();
	unsigned n_layers = LEVELS.size();

	vector<bool> best_path;
	vector<bool> local_path;
	best_path.resize(n_layers, false);
	local_path.resize(n_layers, false);

	unsigned actual_move_in_path = n_layers;
	double actual_cost_in_path = INTMAX_MAX;
	unsigned best_move_in_path = n_layers;
	double best_cost_in_path = INTMAX_MAX;

	S->copy_HDAG(source);
	S->set_cost(best.get_cost());

	for (unsigned iter = 0; iter < n_layers; ++iter) {

		actual_cost_in_path = INTMAX_MAX;
		actual_move_in_path = n_layers;
		best_move_in_path = n_layers;

		for (unsigned l = 0; l < n_layers; ++l) {

			if (local_path[l] == false) {

				// FLIP
				actual_cost = S->get_cost();
				double cost_before_flip = actual_cost;

				if (l > 0)
					actual_cost -= compute_cost_to_subtract(*S, l - 1);
				actual_cost -= compute_cost_to_subtract(*S, l);

				S->copy_level(target, l);

				if (l > 0)
					actual_cost += compute_cost_to_add(*S, l - 1);
				actual_cost += compute_cost_to_add(*S, l);

				S->set_cost(actual_cost);

				if (actual_cost < actual_cost_in_path) {
					actual_cost_in_path = actual_cost;
					actual_move_in_path = l;
				}

				if (actual_cost < best_cost_in_path) {

					best_cost_in_path = actual_cost;
					best_move_in_path = l;

				}

				// UNFLIP
				S->copy_level(source, l);
				S->set_cost(cost_before_flip);

			}
		}

		if (actual_move_in_path != n_layers) {
			S->copy_level(target, actual_move_in_path);
			S->set_cost(actual_cost_in_path);
			local_path[actual_move_in_path] = true;

//			cout << " movement = " << actual_move_in_path << endl;

		}
		if (best_move_in_path != n_layers) {
			best_path[best_move_in_path] = true;
		}

	}

	cout << endl;
	unsigned false_remaining = 0;
	for (unsigned l = 0; l < local_path.size(); ++l)
		if (local_path[l] == false)
			false_remaining++;
	cout << " number FALSE = " << false_remaining << endl;
	cout << " TARGET cost = " << target.get_cost() << endl;
	cout << " BUILT  cost = " << S->get_cost() << endl;

	if (best_cost_in_path < best.get_cost()) {
		for (unsigned l = 0; l < n_layers; ++l) {
			if (best_path[l] == true) {
				best.copy_level(target, l);
			}
		}
		best.set_cost(best_cost_in_path);
	}*/

/*	double actual_cost;
	HDAG *S = new HDAG;

	S->copy_HDAG(source);
	S->set_cost(source.get_cost());

	const auto &L_S = S->getLEVELS();
	for (unsigned l = 0; l < L_S.size(); ++l) {
		actual_cost = S->get_cost();
		// Compute cost to subtract for level l-1
		if (l > 0)
			actual_cost -= compute_cost_to_subtract(*S, l - 1);
		// Compute cost to subtract for level l
		actual_cost -= compute_cost_to_subtract(*S, l);
		S->copy_level(target, l);

		if (l == 0)
		S->printHDAG(*S, "HDAG-tabu+PR-level-0-copied.tex");

		// Compute cost to add for level l-1
		if (l > 0)
			actual_cost += compute_cost_to_add(*S, l - 1);
		// Compute cost to add for level l
		actual_cost += compute_cost_to_add(*S, l);
		S->set_cost(actual_cost);
		if (S->get_cost() < best.get_cost()) {
			if (improved) {
				for (unsigned l_ = last_l_improved + 1;
						l_ <= l; ++l_) {
					best.copy_level(*S, l_);
				}
				last_l_improved = l;
				best.set_cost(S->get_cost());
			}
			else {
				for (unsigned l_ = 0; l_ < l; ++l_) {
					best.copy_level(*S, l_);
				}
				last_l_improved = l;
				improved = true;
				best.set_cost(S->get_cost());
			}
		}

	}*/


	// Maintain the source solution in a dummy structure
	HDAG *S = new HDAG;
	S->copy_HDAG(source);
	S->set_cost(source.get_cost());
	double actual_cost = S->get_cost();
	double best_cost = target.get_cost();
	bool improved = false;
	SolutionIMLCM s;

	const auto &L = S->getLEVELS();
	vector<char> best_configuration;
	best_configuration.resize(L.size(), 'S');

//	auto &LEVELS_target = target.getLEVELS();
//	auto &IDs_target = target.getIDs();

//	auto &IDs_source = source.getIDs();
//	auto &Pos_source = source.getPos();

//	cout << endl;

	for (unsigned l = 0; l < L.size(); ++l) {

//		if (l > 0) subtract_cost_in_level(*S, l - 1, actual_cost);
//		subtract_cost_in_level(*S, l, actual_cost);

/*		for (unsigned i = 0; i < L[l].size(); ++i) {

			unsigned id_v_in_target = IDs_target[l][i];
			unsigned pos_v_in_target = i;

			unsigned id_v_in_source = id_v_in_target;
			unsigned pos_v_in_source = Pos_source[l][id_v_in_source];

			if (pos_v_in_target != pos_v_in_source) {

				s.getCost_beforeOneSwap(I, *S, l, pos_v_in_source, pos_v_in_target, actual_cost);
				S->swapPositions(pos_v_in_source, pos_v_in_target, l);
				s.getCost_afterOneSwap(I, *S, l, pos_v_in_source, pos_v_in_target, actual_cost);

			}

		}*/



//		cout << " --> " << actual_cost << endl;

		S->copy_level(target, l);
		actual_cost = s.getCost(*S);


//		if (l > 0) add_cost_in_level(*S, l - 1, actual_cost);
//		add_cost_in_level(*S, l, actual_cost);

		if (actual_cost < best_cost) {

			improved = true;
			best_cost = actual_cost;
			best_configuration[l] = 'T';

		}

	}

	if (improved) {

		for (unsigned l = 0; l < L.size(); ++l) {

			if (best_configuration[l] == 'S') {
				best.copy_level(source, l);
			}
			else if (best_configuration[l] == 'T') {
				best.copy_level(target, l);
			}

		}

		best.set_cost(best_cost);

	}

	return;

}

Path_Relinking::~Path_Relinking() {
	// TODO Auto-generated destructor stub
}

