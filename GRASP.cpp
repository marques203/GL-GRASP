/*
 * GRASP.cpp
 *
 *  Created on: 23 mag 2017
 *      Author: antonio
 */

#include "GRASP.h"
#include <unistd.h>

GRASP::GRASP() :
		k(0), rhomax(0), maxdeg(0), /*minimum_deg_partial(0), maximum_deg_partial(
				0), tau(0),*/ alpha(0), seed(0), max_iter(100), tenure_vertex_position(1),
				ls(), total_in_RCL(0), total_to_assign(
				0), out_f_best(NULL), out_f_complete(NULL), best_cost(
				INT64_MAX), time_to_best(0), ES(nullptr), PR(nullptr) {
	// TODO Auto-generated constructor stub

}

GRASP::~GRASP() {
	// TODO Auto-generated destructor stub
}

void GRASP::localsearch(HDAG &I, HDAG &S) {

	SolutionIMLCM s;
	s.initialize_costsForLevels(I.getLEVELS().size());

	double best_cost = s.getCost(S);
	S.set_cost(best_cost);

	// try to push up an incremental node until the solution improves;
	// otherwise try to push it down until the solution improves, but
	// this second step must guarantee the feasibility position constraint of the original node
	move(I, S, s, best_cost);

	S.set_cost(best_cost);

	set_end_iteration(chrono::system_clock::now());

	int elapsed =
			chrono::duration_cast<chrono::milliseconds>(end_iteration -
					start_algorithm).count();
	double current = (double)((double)elapsed / (double)1000);

	if (best_cost < get_best_cost()) {

		best_solution.set_cost(best_cost);
		set_best_cost(best_cost);
		set_time_to_best(current);
		for (unsigned i = 0; i < S.getLevNumber(); ++i) {
			best_solution.copy_level(S, i);
		}

	}

//	fprintf(out_f_best, "%d %d\n", (int)constr_cost, (int)best_cost);
//	fprintf(out_f_complete, "%d %d\n", (int)constr_cost, (int)best_cost);

	print_complete_data_out_file(best_cost, current);

	delete_all_iterative_structure(I);

	return;

}

void GRASP::allocate_and_build_all_structure(HDAG &I) {

	auto LevelNumber = I.getLevNumber();

	CL.resize(LevelNumber);
	toAssign.resize(LevelNumber, 0);
	RCL.resize(LevelNumber);
	RhoCL.resize(LevelNumber);
	residual_shift.resize(LevelNumber);
	upmost_inserted_node.resize(LevelNumber);
	OD.resize(LevelNumber);

	for (unsigned l = 0; l < LevelNumber; ++l) {
		unsigned pos_in_level = I.getLevel(l).size();
		CL[l].resize(pos_in_level, 1);
		OD[l].resize(pos_in_level, 0);
		RhoCL[l].resize(pos_in_level, 0);
	}

	return;

}

void GRASP::initialize_all_iteration_structure(HDAG &I) {

	auto LevelNumber = I.getLevNumber();

	setTotalInRCL(0);
	setTotalToAssign(0);

	for (unsigned l = 0; l < LevelNumber; ++l) {
		unsigned pos_in_level = I.getLevel(l).size();
		toAssign[l] = 0;
		for (unsigned i = 0; i < pos_in_level; ++i) {
			CL[l][i] = 1;
			RhoCL[l][i] = 0;
		}
	}

	maxdeg=0;

	return;

}

void GRASP::delete_all_iterative_structure(HDAG &I) {

	auto LevelNumber = I.getLevNumber();
	if (total_in_RCL > 0) {
		for (unsigned l = 0; l < LevelNumber; ++l) {
			while (RCL[l].size() != 0) {
				RCL[l].pop_back();
			}
		}
	}
	if (v_maxoutdeg.size() > 0) {
		while (v_maxoutdeg.size() > 0) {
			v_maxoutdeg.pop_back();
		}
	}

}

void GRASP::build_maxdegree_node_structure(HDAG &I, HDAG &S) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned in_deg = LEVELS[l][i].size();
			unsigned out_deg = B_LEVELS[l][i].size();
			unsigned tot_deg = in_deg + out_deg;
			if (tot_deg > maxdeg)
				maxdeg = tot_deg;
		}
	}


	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned in_deg = LEVELS[l][i].size();
			unsigned out_deg = B_LEVELS[l][i].size();
			unsigned tot_deg = in_deg + out_deg;

			if (tot_deg == maxdeg) {

				pair<unsigned, unsigned> p;
				p.first = l;
				p.second = i;
				v_maxoutdeg.push_back(p);

			}
		}
	}


	return;

}

void GRASP::compute_degree_matrix(HDAG &I) {

	unsigned ln = I.getLevNumber();

	for (unsigned l = 0; l < ln; ++l) {
		unsigned nn = I.getLevel(l).size();
		for (unsigned i = 0; i < nn; ++i) {
			OD[l][i] = I.compute_degree(l, i);
		}
	}
}

pair<unsigned, unsigned> GRASP::extract_maxdegree_node(HDAG &instance, MTRand &rng) {

	std::pair<unsigned, unsigned> v;

	unsigned node = rng.randInt(v_maxoutdeg.size() - 1);

	// Select the element v in v_max_degree and put it
	// in a random position in the solution
	unsigned l_v = v_maxoutdeg[node].first;
	unsigned i_v = v_maxoutdeg[node].second;

	v_maxoutdeg.erase(v_maxoutdeg.begin() + node);

	v.first = l_v;
	v.second = i_v;

	return v;

}

void GRASP::add_link(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) {

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();

	auto &CL = getCL();

	unsigned id_u = IDs[l][i];

	/********************** new feature ***********************/
//	bool possible_minimum_deg_partial_changed = false;
	/**********************************************************/


	// Analyze BS(id_u) and add connections (id_v, id_u) such that id_v \in CL
	if (l > 0) {

		for (unsigned i_ = 0; i_ < LEVELS[l - 1].size(); ++i_) {
			if (CL[l - 1][i_] == 0) {
				for (unsigned j = 0; j < LEVELS[l - 1][i_].size(); ++j) {
					unsigned id_u_ = LEVELS[l - 1][i_][j];
					if (id_u_ == id_u) {
						unsigned id_v = IDs[l - 1][i_];
						unsigned pos_v_in_sol = Pos_S[l - 1][id_v];
						LEVELS_S[l - 1][pos_v_in_sol].push_back(id_u);

						/********************** new feature ***********************/
/*						increase_value_deg_partial(l-1, pos_v_in_sol, 1);
						if (deg_partial[l-1][pos_v_in_sol] > maximum_deg_partial) {
							maximum_deg_partial = deg_partial[l-1][pos_v_in_sol];
						}
						if (minimum_deg_partial ==
								deg_partial[l-1][pos_v_in_sol] - 1) {
							possible_minimum_deg_partial_changed = true;
						}*/
						/**********************************************************/

						break;
					}
				}
			}
		}
	}

	// Analyze FS(id_u) and add connections (id_u, id_v) such that id_v \in CL
	if (l < LEVELS.size() - 1) {

		for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
			unsigned id_v = LEVELS[l][i][j];
			unsigned pos_v = Pos[l + 1][id_v];
			if (CL[l + 1][pos_v] == 0) {
				LEVELS_S[l][bc].push_back(id_v);

				/********************** new feature ***********************/
/*				increase_value_deg_partial(l, bc, 1);
				if (deg_partial[l][bc] > maximum_deg_partial) {
					maximum_deg_partial = deg_partial[l][bc];
				}*/
				/**********************************************************/
			}
		}
	}

	/********************** new feature ***********************/
/*	if (possible_minimum_deg_partial_changed) {
		for (unsigned l = 0; l < LEVELS_S.size(); ++l) {
			for (unsigned i = 0; i < LEVELS_S[l].size(); ++i) {
				if (CL[l][i] == 0) {
					if (minimum_deg_partial >
					deg_partial[l][i])
						minimum_deg_partial = deg_partial[l][i];
				}
			}
		}
	}*/
	/**********************************************************/

	return;

}

void GRASP::update_rho(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	if (l < LEVELS.size() - 1) {
		for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
			unsigned id_ = LEVELS[l][i][j];
			unsigned pos_id_ = Pos[l + 1][id_];
			if (CL[l + 1][pos_id_] == 1) {
				RhoCL[l + 1][pos_id_]++;

				if (RhoCL[l + 1][pos_id_] > rhomax)
					rhomax = RhoCL[l + 1][pos_id_];

				if (!is_already_in_RCL(l + 1, pos_id_)) {

					/********************** new feature ***********************/
					if (RhoCL[l + 1][pos_id_] >= floor(alpha * rhomax)) {
					//if (RhoCL[l+1][pos_id_] >= get_tau()) {
					/**********************************************************/

						pair<unsigned, unsigned> tmp;
						tmp.first = l + 1;
						tmp.second = pos_id_;
						RCL[l + 1].push_back(tmp);
						total_in_RCL++;
					}
				}
			}
		}
	}
	if (l > 0) {
		for (unsigned j = 0; j < B_LEVELS[l][i].size(); ++j) {
			unsigned id_ = B_LEVELS[l][i][j];
			unsigned pos_id_ = Pos[l - 1][id_];
			if (CL[l - 1][pos_id_] == 1) {
				RhoCL[l - 1][pos_id_]++;

				if (RhoCL[l - 1][pos_id_] > rhomax)
					rhomax = RhoCL[l - 1][pos_id_];

				if (!is_already_in_RCL(l - 1, pos_id_)) {

					/********************** new feature ***********************/
					if (RhoCL[l - 1][pos_id_] >= floor(alpha * rhomax)) {
					//if (RhoCL[l-1][pos_id_] >= get_tau()) {
					/**********************************************************/

						pair<unsigned, unsigned> tmp;
						tmp.first = l - 1;
						tmp.second = pos_id_;
						RCL[l - 1].push_back(tmp);
						total_in_RCL++;

					}
				}
			}
		}
	}

	return;

}

void GRASP::update_rhomax(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	if (l > 0) {
		for (unsigned j = 0; j < B_LEVELS[l][i].size(); ++j) {
			unsigned id_ = B_LEVELS[l][i][j];
			unsigned pos_id_ = Pos[l - 1][id_];
			if (CL[l - 1][pos_id_] == 1) {
				if (rhomax < RhoCL[l - 1][pos_id_] + 1) {
					rhomax = RhoCL[l - 1][pos_id_] + 1;
				}
			}
		}
	}

	if (l < LEVELS.size() - 1) {
		for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
			unsigned id_ = LEVELS[l][i][j];
			unsigned pos_id_ = Pos[l + 1][id_];
			if (CL[l + 1][pos_id_] == 1) {
				if (rhomax < RhoCL[l + 1][pos_id_] + 1) {
					rhomax = RhoCL[l + 1][pos_id_] + 1;
				}
			}
		}
	}

	return;

}

void GRASP::update_RCL(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &LEVELS = I.getLEVELS();
	unsigned new_rhomax = 0;

	// rhomax could be changed
	if (RhoCL[l][i] == rhomax) {
		for (unsigned l = 0; l < LEVELS.size(); ++l) {
			for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
				if (CL[l][i] == 1) {
					if (RhoCL[l][i] > new_rhomax)
						new_rhomax = RhoCL[l][i];
				}
			}
		}
	}

	rhomax = new_rhomax;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			if (CL[l][i] == 1) {
				if (!is_already_in_RCL(l, i)) {
					if (RhoCL[l][i] >= floor(alpha * rhomax)) {
						pair<unsigned, unsigned> tmp;
						tmp.first = l;
						tmp.second = i;
						RCL[l].push_back(tmp);
						total_in_RCL++;
					}
				}
				else {
					if (RhoCL[l][i] < floor(alpha * rhomax)) {
						RCL[l].erase(RCL[l].begin() + i);
						total_in_RCL--;
					}
				}
			}
		}
	}

}

bool GRASP::is_already_in_RCL(unsigned l, unsigned i) {

	for (unsigned i_ = 0; i_ < RCL[l].size(); ++i_)
		if (RCL[l][i_].second == i)
			return true;

	return false;

}

void GRASP::build_RCL(HDAG &instance, HDAG &solution) {

	auto &CL = getCL();
	auto &RCL = getRCL();

	for (unsigned l = 0; l < CL.size(); ++l) {

		for (unsigned i = 0; i < CL[l].size(); ++i) {

			if (CL[l][i] == 1) {

				unsigned rho = compute_potential_degree(instance, solution, l,
						i);

				RhoCL[l][i] = rho;

				if (rho >= floor(getalpha() * rhomax)) {

					pair<unsigned, unsigned> tmp;
					tmp.first = l;
					tmp.second = i;

					RCL[l].push_back(tmp);
					getTotalInRCL()++;

				}

			}

		}

	}

}

unsigned GRASP::compute_potential_degree(HDAG &instance, HDAG &solution,
		unsigned l, unsigned i) {

	auto &CL = getCL();

	const auto &LEVELS = instance.getLEVELS();
	const auto &IDs = instance.getIDs();

	unsigned rho = 0;

	// Explore the backward star
	if (l > 0) {

		for (unsigned i_ = 0; i_ < CL[l - 1].size(); ++i_) {

			// Then the element in position [l-1][i_] in instance
			// is in solution
			if (CL[l - 1][i_] == 0) {

				// Now explore its Forward Start
				for (unsigned j_ = 0; j_ < LEVELS[l - 1][i_].size(); ++j_) {

					// Check if the considered node is in its Forward Star
					if (LEVELS[l - 1][i_][j_] == IDs[l][i]) {

						rho++;
						break;

					}

				}

			}

		}

	}

	// Explore the forward star
	if (l < CL.size() - 1) {

		// Scan the forward star of the considered node
		for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

			// For all the elements of the considered node
			// that are in solution, increment the value of rho
			if (is_in_CL(instance, l + 1, LEVELS[l][i][j])) {

				rho++;

			}
		}

	}

	return rho;

}

bool GRASP::is_in_CL(HDAG &instance, unsigned l, unsigned id) {

	unsigned pos_in_I = instance.getPos()[l][id];
	if (getCL()[l][pos_in_I] == 0)
		return true;

	return false;
}

pair<unsigned, unsigned> GRASP::extract_from_RCL(MTRand &rng) {

	auto &RCL = getRCL();

	pair<unsigned, unsigned> tmp;

	unsigned l, i;

	while (1) {

		l = rng.randInt(RCL.size() - 1);
		if (RCL[l].size() > 0)
			break;

	}

	i = rng.randInt(RCL[l].size() - 1);
	tmp.first = RCL[l][i].first;
	tmp.second = RCL[l][i].second;
	RCL[l].erase(RCL[l].begin() + i);

	getTotalInRCL()--;

	return tmp	;

}

bool GRASP::is_position_free(HDAG &solution, unsigned l, unsigned bc) {

	const auto &sol_IDs = solution.getIDs();
	const auto &sol_LEVELS = solution.getLEVELS();

	unsigned vn = sol_LEVELS[l].size();

	if (sol_IDs[l][bc] != vn)
		return false;

	return true;

}

void GRASP::compute_rhomax(HDAG &I, HDAG &S) {

	rhomax = 0;

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	for (unsigned l = 1; l < LEVELS.size() - 1; ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			if (CL[l][i] == 1) {
				unsigned rho = 0;
				for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
					unsigned id_ = LEVELS[l][i][j];
					unsigned pos_id_ = Pos[l + 1][id_];
					if (CL[l + 1][pos_id_] == 0) {
						rho++;
					}
				}
				for (unsigned i_ = 0; i_ < B_LEVELS[l][i].size(); ++i_) {
					unsigned id_ = B_LEVELS[l][i][i_];
					unsigned pos_id_ = Pos[l - 1][id_];
					if (CL[l - 1][pos_id_] == 0) {
						rho++;
					}
				}
				rhomax = max(rho, rhomax);
			}
		}
	}

	// Only for the level-0
	for (unsigned i = 0; i < LEVELS[0].size(); ++i) {
		if (CL[0][i] == 1) {
			unsigned rho = 0;
			for (unsigned j = 0; j < LEVELS[0][i].size(); ++j) {
				unsigned id_ = LEVELS[0][i][j];
				unsigned pos_id_ = Pos[1][id_];
				if (CL[1][pos_id_] == 0) {
					rho++;
				}
			}
			rhomax = max(rho, rhomax);
		}
	}

	// Only for the last level
	unsigned last_l = LEVELS.size() - 1;
	for (unsigned i = 0; i < LEVELS[last_l].size(); ++i) {
		if (CL[last_l][i] == 1) {
			unsigned rho = 0;
			for (unsigned i_ = 0; i_ < B_LEVELS[last_l][i].size(); ++i_) {
				unsigned id_ = B_LEVELS[last_l][i][i_];
				unsigned pos_id_ = Pos[last_l - 1][id_];
				if (CL[last_l - 1][pos_id_] == 0) {
					rho++;
				}
			}
			rhomax = max(rho, rhomax);
		}
	}

	return;

}

void GRASP::phase_1(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost) {

	bool improve = true;
	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();
	double cost = best_cost;
	double prev_cost = cost;
	auto &INCREMENTAL = I.getINCREMENTAL();

	while (improve) {
		improve = false;
		for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
			cost = best_cost;
			for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
				unsigned pos = Pos_S[l][INCREMENTAL[l][i]];
				unsigned best_flip = LEVELS_S[l].size();
				for (unsigned i_ = 0; i_ < INCREMENTAL[l].size(); ++i_) {
					if (INCREMENTAL[l][i] != INCREMENTAL[l][i_]) {
						unsigned pos_ = Pos_S[l][INCREMENTAL[l][i_]];
						prev_cost = cost;

						cost = s.getCost_beforeOneSwap(I, S, l, pos, pos_, cost);
						S.swapPositions(pos, pos_, l);
						cost = s.getCost_afterOneSwap(I, S, l, pos, pos_, cost);

						if (cost < best_cost) {
							best_cost = cost;
							best_flip = pos_;
							if (l > 0)
								s.get_bestsCostsForLevels()[l - 1] =
										s.get_costsForLevels()[l - 1];
							s.get_bestsCostsForLevels()[l] =
									s.get_costsForLevels()[l];
						}
						S.swapPositions(pos_, pos, l);
						cost = prev_cost;
					}
				}
				if (best_flip != LEVELS_S[l].size()) {
					S.swapPositions(pos, best_flip, l);
					improve = true;
				}
				best_flip = LEVELS_S[l].size();
				cost = best_cost;
			}
		}
	}

	return;
}

void GRASP::move(HDAG &I, HDAG &S, 	SolutionIMLCM &s, double &best_cost) {

//	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();

	auto &INCREMENTAL = I.getINCREMENTAL();

	if (strcmp(getlocalsearch(), "best") == 0) {
		// best improvement local search
		phase_1(I, S, s, best_cost);
		unsigned best_shift_up;
		unsigned best_shift_dw;
		bool improve = true;
		while (improve) {
			double old_cost = best_cost;
			for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
				for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
					best_shift_up = find_best_shift_up(I, S,
							Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
							improve);
					shift_down(S, Pos_S[l][INCREMENTAL[l][i]], best_shift_up, l);

					if (old_cost == best_cost) {
						best_shift_dw = find_best_shift_down(I, S,
								Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
								improve);
						shift_up(S, Pos_S[l][INCREMENTAL[l][i]], best_shift_dw,
								l);
					}
				}
			}
			if (old_cost == best_cost) {
				improve = false;
			}
		}
	} else {
		// first improvement local search
		unsigned first_shift_up;
		unsigned first_shift_dw;
		bool improve = true;
		if (strcmp(ls, "first-prepro") == 0)
			phase_1(I, S, s, best_cost);
		unsigned number_iteration=0;
		while (improve) {
			if (number_iteration >= 100) break;
			double old_cost = best_cost;

			for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
				for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
					first_shift_up = find_first_shift_up(I, S,
							Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
							improve);
					shift_down(S, Pos_S[l][INCREMENTAL[l][i]], first_shift_up, l);
					if (old_cost == best_cost) {
						first_shift_dw = find_first_shift_down(I, S,
								Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
								improve);
						shift_up(S, Pos_S[l][INCREMENTAL[l][i]], first_shift_dw,
								l);
					}
				}
			}
			if (old_cost == best_cost) {
				improve = false;
			}

			number_iteration++;

		}
	}

	return;

}

unsigned GRASP::find_best_shift_up(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve) {

	auto &LEVELS_S = S.getLEVELS();
	unsigned best_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {
		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);
//		cost = s.getCost_afterSwap(S, l, cost);
		if (cost < best_cost) {
			best_cost = cost;
			best_shift = i + 1;
			if (l > 0)
				s.get_bestsCostsForLevels()[l - 1] =
			    s.get_costsForLevels()[l - 1];
			s.get_bestsCostsForLevels()[l] =
			s.get_costsForLevels()[l];
			improve = true;
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

unsigned GRASP::find_first_shift_up(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve) {

	auto &LEVELS_S = S.getLEVELS();
	unsigned first_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {
		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);
		//cost = s.getCost_afterSwap(S, l, cost);
		if (cost < best_cost) {
			best_cost = cost;
			first_shift = i + 1;
			if (l > 0)
				s.get_bestsCostsForLevels()[l - 1] =
			    s.get_costsForLevels()[l - 1];
			s.get_bestsCostsForLevels()[l] =
			s.get_costsForLevels()[l];
			return first_shift;
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return first_shift;

}


void GRASP::shift_down(HDAG &S, unsigned up_limit, unsigned dw_limit, unsigned l) {

	for (unsigned i = up_limit; i > dw_limit; --i)
		S.swapPositions(i, i-1, l);

}

unsigned GRASP::find_best_shift_down(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	unsigned best_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {
				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
				//cost = s.getCost_afterSwap(S, l, cost);
				if (cost < best_cost) {
					best_cost = cost;
					best_shift = i - 1;
					if (l > 0)
						s.get_bestsCostsForLevels()[l - 1] =
					    s.get_costsForLevels()[l - 1];
					s.get_bestsCostsForLevels()[l] =
					s.get_costsForLevels()[l];
					improve = true;
				}

			}
			else {
				if (l > 0)
					s.get_costsForLevels()[l - 1] =
				    s.get_bestsCostsForLevels()[l - 1];
				s.get_costsForLevels()[l] =
				s.get_bestsCostsForLevels()[l];
				return best_shift;

			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
			//cost = s.getCost_afterSwap(S, l, cost);
			if (cost < best_cost) {
				best_cost = cost;
				best_shift = i - 1;
				if (l > 0)
					s.get_bestsCostsForLevels()[l - 1] =
				    s.get_costsForLevels()[l - 1];
				s.get_bestsCostsForLevels()[l] =
				s.get_costsForLevels()[l];
				improve = true;
			}
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

unsigned GRASP::find_first_shift_down(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	double cost = best_cost;

	for (unsigned i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {
				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
				//cost = s.getCost_afterSwap(S, l, cost);
				if (cost < best_cost) {
					best_cost = cost;
					if (l > 0)
						s.get_bestsCostsForLevels()[l - 1] =
					    s.get_costsForLevels()[l - 1];
					s.get_bestsCostsForLevels()[l] =
					s.get_costsForLevels()[l];
					return i - 1;
				}

			}
			else {
				if (l > 0)
					s.get_costsForLevels()[l - 1] =
				    s.get_bestsCostsForLevels()[l - 1];
				s.get_costsForLevels()[l] =
				s.get_bestsCostsForLevels()[l];

				return init;
			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
			//cost = s.getCost_afterSwap(S, l, cost);
			if (cost < best_cost) {
				best_cost = cost;
				if (l > 0)
					s.get_bestsCostsForLevels()[l - 1] =
				    s.get_costsForLevels()[l - 1];
				s.get_bestsCostsForLevels()[l] =
				s.get_costsForLevels()[l];
				return i - 1;
			}
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return init;

}

void GRASP::shift_up(HDAG &S, unsigned dw_limit, unsigned up_limit, unsigned l) {

	for (unsigned i = dw_limit; i < up_limit; ++i)
		S.swapPositions(i, i + 1, l);

}

void GRASP::allocate_best_solution(HDAG &I) {

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &best_solution_IDs = best_solution.getIDs();
	auto &best_solution_Pos = best_solution.getPos();
	auto &best_solution_Os = best_solution.getOs();

	unsigned ln, vn;

	ln = I.getLevNumber();

	best_solution.allocateLEVELS(ln);
	best_solution.allocateIDs(ln);
	best_solution.allocateOs(ln);
	best_solution.allocatePos(ln);

	for (unsigned l = 0; l < ln; ++l) {

		vn = LEVELS[l].size();

		best_solution.allocateLevel(l, vn);
		best_solution.allocateIDs(l, vn);
		best_solution.allocateOs(l, vn);
		best_solution.allocatePos(l, vn);

		for (unsigned i = 0; i < vn; ++i) {

			unsigned id_u = IDs[l][i];
			best_solution_IDs[l][i] = id_u;
			best_solution_Pos[l][id_u] = i;
			best_solution_Os[l][i] = Os[l][i];

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				best_solution.getLEVELS()[l][i].push_back(LEVELS[l][i][j]);

			}

		}

	}

	return;

}

void GRASP::copy_best_solution(HDAG &S) {

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	auto &best_solution_IDs = best_solution.getIDs();
	auto &best_solution_Pos = best_solution.getPos();

	for (unsigned l = 0; l < LEVELS_S.size(); ++l) {
		for (unsigned i = 0; i < LEVELS_S[l].size(); ++i) {
			if (best_solution_IDs[l][i] != IDs_S[l][i]) {
				unsigned id_S = IDs_S[l][i];
				unsigned Pos_bs = best_solution_Pos[l][id_S];
				best_solution.swapPositions(i, Pos_bs, l);
			}
		}
	}

	return;

}

/**************************************** FOR TABU ***************************************/
void GRASP::build_TABU_POSITION_MATRIX(HDAG &I) {
	TABU_POSITION_MATRIX.resize(I.getLEVELS().size());
	for (unsigned l = 0; l < I.getLEVELS().size(); ++l) {
		TABU_POSITION_MATRIX[l].resize(I.getLEVELS()[l].size());
		for (unsigned i = 0; i < I.getLEVELS()[l].size(); ++i) {
			TABU_POSITION_MATRIX[l][i].resize(I.getLEVELS()[l].size(), 0);
		}
	}
}

void GRASP::build_TABU_VERTEX_MATRIX(HDAG &I) {
	TABU_VERTEX_MATRIX.resize(I.getLEVELS().size());
	for (unsigned l = 0; l < I.getLEVELS().size(); ++l) {
		TABU_VERTEX_MATRIX[l].resize(I.getLEVELS()[l].size(), 0);
	}
}


void GRASP::update_TABU_POSITION_MATRIX(HDAG &I, HDAG &S, unsigned iter) {

	auto &LEVELS_S = S.getLEVELS();

	auto &IDs_S = S.getIDs();

	auto &Pos = I.getPos();

	auto &Os_S = S.getOs();

	for (unsigned l = 0; l < LEVELS_S.size(); ++l) {
		for (unsigned i = 0; i < LEVELS_S[l].size(); ++i) {
			if (Os_S[l][i] == 0) {
				unsigned id = IDs_S[l][i];
				unsigned original_pos = Pos[l][id];

				TABU_POSITION_MATRIX[l][original_pos][i] = iter;

			}
		}
	}

	return;

}

void GRASP::phase_1_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost, unsigned iter) {

	bool improve = true;
	auto &Pos = I.getPos();
	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();
	double cost = best_cost;
	double prev_cost = cost;
	auto &INCREMENTAL = I.getINCREMENTAL();

	while (improve) {
		improve = false;
		for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
			cost = best_cost;
			for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
				unsigned pos = Pos_S[l][INCREMENTAL[l][i]];
				unsigned original_pos_i = Pos[l][INCREMENTAL[l][i]];
				unsigned best_flip = LEVELS_S[l].size();
				for (unsigned i_ = 0; i_ < INCREMENTAL[l].size(); ++i_) {
					if (INCREMENTAL[l][i] != INCREMENTAL[l][i_]) {

						unsigned pos_ = Pos_S[l][INCREMENTAL[l][i_]];
						unsigned original_pos_i_ = Pos[l][INCREMENTAL[l][i_]];

						prev_cost = cost;

						cost = s.getCost_beforeOneSwap(I, S, l, pos, pos_, cost);
						S.swapPositions(pos, pos_, l);
						cost = s.getCost_afterOneSwap(I, S, l, pos, pos_, cost);

						// Do the swap only if:
						// the new cost is better and
						// the new positions are not TABU positions for the nodes involved in the swap
						if ( (cost < best_cost) &&
							 !is_TABU_POSITION(I, S, l, original_pos_i, pos_, iter) &&
							 !is_TABU_POSITION(I, S, l, original_pos_i_, pos, iter) ) {

							best_cost = cost;
							best_flip = pos_;
							if (l > 0)
								s.get_bestsCostsForLevels()[l - 1] =
									s.get_costsForLevels()[l - 1];
							s.get_bestsCostsForLevels()[l] =
								s.get_costsForLevels()[l];

						}
						S.swapPositions(pos_, pos, l);
						cost = prev_cost;

					}

				}
				if (best_flip != LEVELS_S[l].size()) {
					S.swapPositions(pos, best_flip, l);
					improve = true;
				}
				best_flip = LEVELS_S[l].size();
				cost = best_cost;
			}
		}
	}

	return;

}

unsigned GRASP::find_best_shift_up_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &Pos = I.getPos();
	auto &IDs_S = S.getIDs();

	auto &LEVELS_S = S.getLEVELS();
	unsigned best_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {

		unsigned id_i_to_shift = IDs_S[l][i];
		unsigned id_i_shifted = IDs_S[l][i + 1];
		unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
		unsigned pos_i_shifted = Pos[l][id_i_shifted];

		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);

		// Consider the shift if only if:
		// the new cost is better and
		// the new positions are not TABU positions for the nodes involved in the swap
		if ( cost < best_cost &&
			 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i + 1, iter) &&
			 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {

			best_cost = cost;
			best_shift = i + 1;
			if (l > 0)
				s.get_bestsCostsForLevels()[l - 1] =
			    s.get_costsForLevels()[l - 1];
			s.get_bestsCostsForLevels()[l] =
			s.get_costsForLevels()[l];
			improve = true;

		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

unsigned GRASP::find_first_shift_up_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &Pos = I.getPos();
	auto &IDs_S = I.getIDs();

	auto &LEVELS_S = S.getLEVELS();
	unsigned first_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {
		unsigned id_i_to_shift = IDs_S[l][i];
		unsigned id_i_shifted = IDs_S[l][i + 1];
		unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
		unsigned pos_i_shifted = Pos[l][id_i_shifted];

		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);

		// Consider the shift if only if:
		// the new cost is better and
		// the new positions are not TABU positions for the nodes involved in the swap
		if ( cost < best_cost &&
			 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i + 1, iter) &&
			 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {
			best_cost = cost;
			first_shift = i + 1;
			if (l > 0)
				s.get_bestsCostsForLevels()[l - 1] =
			    s.get_costsForLevels()[l - 1];
			s.get_bestsCostsForLevels()[l] =
			s.get_costsForLevels()[l];
			return first_shift;
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return first_shift;

}

unsigned GRASP::find_best_shift_down_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	unsigned best_shift = init;
	double cost = best_cost;

	for (unsigned i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {
				unsigned id_i_to_shift = IDs_S[l][i];
				unsigned id_i_shifted = IDs_S[l][i - 1];
				unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
				unsigned pos_i_shifted = Pos[l][id_i_shifted];

				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);

				// Consider the shift if only if:
				// the new cost is better and
				// the new positions are not TABU positions for the nodes involved in the swap
				if ( cost < best_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {
					best_cost = cost;
					best_shift = i - 1;
					if (l > 0)
						s.get_bestsCostsForLevels()[l - 1] =
					    s.get_costsForLevels()[l - 1];
					s.get_bestsCostsForLevels()[l] =
					s.get_costsForLevels()[l];
					improve = true;
				}

			}
			else {
				if (l > 0)
					s.get_costsForLevels()[l - 1] =
				    s.get_bestsCostsForLevels()[l - 1];
				s.get_costsForLevels()[l] =
				s.get_bestsCostsForLevels()[l];
				return best_shift;

			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
			//cost = s.getCost_afterSwap(S, l, cost);
			unsigned id_i_to_shift = IDs_S[l][i];
			unsigned id_i_shifted = IDs_S[l][i - 1];
			unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
			unsigned pos_i_shifted = Pos[l][id_i_shifted];

			if ( cost < best_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter) ) {
				best_cost = cost;

				best_shift = i - 1;
				if (l > 0)
					s.get_bestsCostsForLevels()[l - 1] =
				    s.get_costsForLevels()[l - 1];
				s.get_bestsCostsForLevels()[l] =
				s.get_costsForLevels()[l];
				improve = true;
			}
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

unsigned GRASP::find_first_shift_down_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	double cost = best_cost;

	for (unsigned i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {

				unsigned id_i_to_shift = IDs_S[l][i];
				unsigned id_i_shifted = IDs_S[l][i - 1];
				unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
				unsigned pos_i_shifted = Pos[l][id_i_shifted];

				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);

				// Consider the shift if only if:
				// the new cost is better and
				// the new positions are not TABU positions for the nodes involved in the swap
				if ( cost < best_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {
					best_cost = cost;
					if (l > 0)
						s.get_bestsCostsForLevels()[l - 1] =
					    s.get_costsForLevels()[l - 1];
					s.get_bestsCostsForLevels()[l] =
					s.get_costsForLevels()[l];
					return i - 1;
				}

			}
			else {
				if (l > 0)
					s.get_costsForLevels()[l - 1] =
				    s.get_bestsCostsForLevels()[l - 1];
				s.get_costsForLevels()[l] =
				s.get_bestsCostsForLevels()[l];

				return init;
			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
			//cost = s.getCost_afterSwap(S, l, cost);

			unsigned id_i_to_shift = IDs_S[l][i];
			unsigned id_i_shifted = IDs_S[l][i - 1];
			unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
			unsigned pos_i_shifted = Pos[l][id_i_shifted];

			if ( cost < best_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter) ) {

				best_cost = cost;
				if (l > 0)
					s.get_bestsCostsForLevels()[l - 1] =
				    s.get_costsForLevels()[l - 1];
				s.get_bestsCostsForLevels()[l] =
				s.get_costsForLevels()[l];
				return i - 1;
			}
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return init;

}

void GRASP::move_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost,
		unsigned iter) {

	//	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();

	auto &INCREMENTAL = I.getINCREMENTAL();

	if (strcmp(getlocalsearch(), "best") == 0) {

		// best improvement local search
		phase_1_TABU(I, S, s, best_cost, iter);
		unsigned best_shift_up;
		unsigned best_shift_dw;
		bool improve = true;
		while (improve) {
			double old_cost = best_cost;
			for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
				for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {

					best_shift_up = find_best_shift_up_TABU(I, S,
							Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
							improve, iter);
					shift_down(S, Pos_S[l][INCREMENTAL[l][i]], best_shift_up,
							l);

					if (old_cost == best_cost) {
						best_shift_dw = find_best_shift_down_TABU(I, S,
								Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
								improve, iter);
						shift_up(S, Pos_S[l][INCREMENTAL[l][i]], best_shift_dw,
								l);

					}
				}
			}
			if (old_cost == best_cost) {
				improve = false;
			}
		}

	} else {
		// first improvement local search
		unsigned first_shift_up;
		unsigned first_shift_dw;
		bool improve = true;
		if (strcmp(ls, "first-prepro") == 0)
			phase_1_TABU(I, S, s, best_cost, iter);
		unsigned number_iteration = 0;
		while (improve) {
			if (number_iteration >= 100)
				break;
			double old_cost = best_cost;

			for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
				for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
					first_shift_up = find_first_shift_up_TABU(I, S,
							Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
							improve, iter);
					shift_down(S, Pos_S[l][INCREMENTAL[l][i]], first_shift_up,
							l);
					if (old_cost == best_cost) {
						first_shift_dw = find_first_shift_down_TABU(I, S,
								Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s,
								improve, iter);
						shift_up(S, Pos_S[l][INCREMENTAL[l][i]], first_shift_dw,
								l);
					}
				}
			}
			if (old_cost == best_cost) {
				improve = false;
			}

			number_iteration++;

		}
	}

	return;

}

unsigned GRASP::find_shift_up_to_restablish_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &Pos = I.getPos();
	auto &IDs_S = S.getIDs();

	auto &LEVELS_S = S.getLEVELS();
	unsigned best_shift = init;
	double initial_cost = INT_MAX;
	double cost = best_cost;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {

		unsigned id_i_to_shift = IDs_S[l][i];
		unsigned id_i_shifted = IDs_S[l][i + 1];
		unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
		unsigned pos_i_shifted = Pos[l][id_i_shifted];

		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);

		// Consider the shift if only if:
		// the new cost is better and
		// the new positions are not TABU positions for the nodes involved in the swap
		if ( cost < initial_cost &&
			 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i + 1, iter) &&
			 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {

			initial_cost = cost;
			best_cost = cost;
			best_shift = i + 1;
			if (l > 0)
				s.get_bestsCostsForLevels()[l - 1] =
			    s.get_costsForLevels()[l - 1];
			s.get_bestsCostsForLevels()[l] =
			s.get_costsForLevels()[l];
			improve = true;

		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

unsigned GRASP::find_shift_down_to_restablish_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	unsigned best_shift = init;
	double initial_cost = INT_MAX;
	double cost = best_cost;

	for (unsigned i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {
				unsigned id_i_to_shift = IDs_S[l][i];
				unsigned id_i_shifted = IDs_S[l][i - 1];
				unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
				unsigned pos_i_shifted = Pos[l][id_i_shifted];

				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);

				// Consider the shift if only if:
				// the new cost is better and
				// the new positions are not TABU positions for the nodes involved in the swap
				if ( cost < initial_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter)) {

					best_cost = cost;
					initial_cost = cost;

					best_shift = i - 1;
					if (l > 0)
						s.get_bestsCostsForLevels()[l - 1] =
					    s.get_costsForLevels()[l - 1];
					s.get_bestsCostsForLevels()[l] =
					s.get_costsForLevels()[l];
					improve = true;
				}

			}
			else {
				if (l > 0)
					s.get_costsForLevels()[l - 1] =
				    s.get_bestsCostsForLevels()[l - 1];
				s.get_costsForLevels()[l] =
				s.get_bestsCostsForLevels()[l];
				return best_shift;

			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);
			//cost = s.getCost_afterSwap(S, l, cost);
			unsigned id_i_to_shift = IDs_S[l][i];
			unsigned id_i_shifted = IDs_S[l][i - 1];
			unsigned pos_i_to_shift = Pos[l][id_i_to_shift];
			unsigned pos_i_shifted = Pos[l][id_i_shifted];

			if ( cost < initial_cost &&
					 !is_TABU_POSITION(I, S, l, pos_i_to_shift, i - 1, iter) &&
					 !is_TABU_POSITION(I, S, l, pos_i_shifted, i, iter) ) {

				initial_cost = cost;
				best_cost = cost;

				best_shift = i - 1;
				if (l > 0)
					s.get_bestsCostsForLevels()[l - 1] =
				    s.get_costsForLevels()[l - 1];
				s.get_bestsCostsForLevels()[l] =
				s.get_costsForLevels()[l];
				improve = true;
			}
		}
	}

	if (l > 0)
		s.get_costsForLevels()[l - 1] =
	    s.get_bestsCostsForLevels()[l - 1];
	s.get_costsForLevels()[l] =
	s.get_bestsCostsForLevels()[l];

	return best_shift;

}

void GRASP::move_to_restablish_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s,
		double &best_cost, unsigned iter) {

	//	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();

	auto &INCREMENTAL = I.getINCREMENTAL();

	// first improvement local search
	unsigned up;
	unsigned dw;
	bool improve = true;
	double old_cost = best_cost;

	for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
		for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {
			up = find_shift_up_to_restablish_TABU(I, S,
					Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s, improve,
					iter);
			shift_down(S, Pos_S[l][INCREMENTAL[l][i]], up, l);
			if (old_cost == best_cost) {
				dw = find_shift_down_to_restablish_TABU(I, S,
						Pos_S[l][INCREMENTAL[l][i]], l, best_cost, s, improve,
						iter);
				shift_up(S, Pos_S[l][INCREMENTAL[l][i]], dw, l);
			}
		}
	}

	return;

}

void GRASP::localsearch_TABU(HDAG &I, HDAG &S, unsigned iter) {

	SolutionIMLCM s;
	s.initialize_costsForLevels(I.getLEVELS().size());

	double actual_best_cost = s.getCost(S);
	S.set_cost(actual_best_cost);

	// try to push up an incremental node until the solution improves;
	// otherwise try to push it down until the solution improves, but
	// this second step must guarantee the feasibility position constraint of the original node
	move_TABU_2(I, S, s, actual_best_cost, iter);

	S.set_cost(actual_best_cost);

	set_end_iteration(chrono::system_clock::now());

	int elapsed =
			chrono::duration_cast<chrono::milliseconds>(end_iteration -
					start_algorithm).count();
	double current = (double)((double)elapsed / (double)1000);

	if (actual_best_cost < get_best_cost()) {

		best_solution.set_cost(actual_best_cost);
		set_best_cost(actual_best_cost);
		set_time_to_best(current);

	}

	print_complete_data_out_file(actual_best_cost, current);

	delete_all_iterative_structure(I);

	update_TABU_POSITION_MATRIX(I, S, iter);

	return;

}



/*********************************** TABU 2 *********************************************/
unsigned GRASP::find_best_shift_up_TABU_2(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &cost_tabu, double &current_cost, SolutionIMLCM &s) {

//	cout << " cost computed before " << s.getCost(S) << endl;

	auto &LEVELS_S = S.getLEVELS();
	auto &Pos_S = S.getPos();

	unsigned best_shift = init;
	double cost = cost_tabu;

	unsigned id_node = S.getIDs()[l][init];

//	cout << " start - current cost " << current_cost << endl;
//	cout << "       - tabu cost    " << cost_tabu << endl;

	for (unsigned i = init; i < LEVELS_S[l].size() - 1; ++i) {

		cost = s.getCost_beforeOneShift(I, S, l, i + 1, i, cost);
		S.swapPositions(i, i + 1, l);
		cost = s.getCost_afterOneShift(I, S, l, i + 1, i, cost);

//		cout << "         dynamic cost " << cost << " shift " << i + 1 << " cost computed during " << s.getCost(S) << endl;

		if (current_cost > cost) {

			current_cost = cost;
			best_shift = i + 1;

		}

		if ( i == LEVELS_S[l].size() - 1 ) break;

	}

	shift_down(S, Pos_S[l][id_node], init, l);

/*	cout << "         init  " << init << endl;
	cout << "         ends  " << Pos_S[l][id_node] << endl;

	cout << "         shift " << best_shift << endl;

	cout << " cost computed after " << s.getCost(S) << endl;

	getchar();*/

	return best_shift;

}

unsigned GRASP::find_best_shift_down_TABU_2(HDAG &I, HDAG &S, unsigned init, unsigned l,
		double &cost_tabu, double &current_cost, SolutionIMLCM &s) {

	auto &IDs_S = S.getIDs();
	auto &Pos = I.getPos();
	auto &Pos_S = S.getPos();

	unsigned id_node = S.getIDs()[l][init];

	unsigned best_shift = init;
	double cost = cost_tabu;

	unsigned i;

	for (i = init; i > 0; --i) {
		if (S.isOriginalNode(i - 1, l)) {
			unsigned id = IDs_S[l][i - 1];
			unsigned origin_p = Pos[l][id];
			unsigned max_feasible_p = origin_p + k;

			if (i <= max_feasible_p) {

				cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
				S.swapPositions(i, i - 1, l);
				cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);

				if (current_cost > cost) {

					current_cost = cost;
					best_shift = i - 1;

				}

			}
			else {

				shift_up(S, Pos_S[l][id_node], init, l);

				return best_shift;

			}
		}
		else {
			cost = s.getCost_beforeOneShift(I, S, l, i, i - 1, cost);
			S.swapPositions(i, i - 1, l);
			cost = s.getCost_afterOneShift(I, S, l, i, i - 1, cost);

			if (current_cost > cost) {

				current_cost = cost;
				best_shift = i - 1;

			}
		}
	}

	shift_up(S, Pos_S[l][id_node], init, l);

	return best_shift;

}

void GRASP::move_TABU_2(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost,
		unsigned iter) {

	phase_1(I, S, s, best_cost);

	auto &Pos = I.getPos();
	auto &Pos_S = S.getPos();

	auto &INCREMENTAL = I.getINCREMENTAL();

	double cost_tabu = best_cost;

	unsigned best_shift_up;
	unsigned best_shift_dw;

	for (unsigned it_tabu = tenure_vertex_position + 1; it_tabu < 100000; it_tabu++) {

		for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {

			double current_cost = INT64_MAX;
			int id_vertex_to_move = -1;
			int old_position_vertex = -1;
			int new_position_vertex = -1;
			int original_position_vertex = -1;
			double initial_cost = cost_tabu;

//			cout << " LEVEL " << l << endl;

			for (unsigned i = 0; i < INCREMENTAL[l].size(); ++i) {

				if (!is_TABU_VERTEX(l, Pos[l][INCREMENTAL[l][i]], it_tabu)) {

					unsigned initial_position = Pos_S[l][INCREMENTAL[l][i]];

/*					cout << " BEFORE " << endl;
					cout << " current cost " << current_cost << endl;
					cout << " cost tabu    " << cost_tabu << endl;
					cout << " id v to move " << id_vertex_to_move << endl;*/

					best_shift_up
					= find_best_shift_up_TABU_2(
						I, S, Pos_S[l][INCREMENTAL[l][i]], l, cost_tabu, current_cost, s);

					if (best_shift_up != initial_position) {
						old_position_vertex = initial_position;
						new_position_vertex = best_shift_up;
						id_vertex_to_move = INCREMENTAL[l][i];
						original_position_vertex = Pos[l][id_vertex_to_move];
					}

					cost_tabu = initial_cost;

					best_shift_dw
					= find_best_shift_down_TABU_2(
						I, S, Pos_S[l][INCREMENTAL[l][i]], l, cost_tabu, current_cost, s);

					if (best_shift_dw != initial_position) {
						old_position_vertex = initial_position;
						new_position_vertex = best_shift_dw;
						id_vertex_to_move = INCREMENTAL[l][i];
						original_position_vertex = Pos[l][id_vertex_to_move];
					}

					cost_tabu = initial_cost;

/*					cout << " AFTER " << endl;
					cout << " current cost " << current_cost << endl;
					cout << " cost tabu    " << cost_tabu << endl;
					cout << " id v to move " << id_vertex_to_move << endl;

					getchar();*/

				}

			}

			if (id_vertex_to_move != -1) {

/*				cout << "    best cost " << best_cost << endl;
				cout << " current cost " << current_cost << endl;
				cout << "   tabu  cost " << cost_tabu << endl;
				getchar();*/


				cost_tabu = current_cost;

				if (old_position_vertex < new_position_vertex) {

					// Shift up
					for (int i = old_position_vertex; i < new_position_vertex; ++i) {
						S.swapPositions(i, i + 1, l);
					}

				}
				else {

					// Shift down
					for (int i = old_position_vertex; i > new_position_vertex; --i) {
						S.swapPositions(i, i - 1, l);
					}

				}

				// Update TABU POSITION MATRIX
				TABU_VERTEX_MATRIX[l][original_position_vertex] = it_tabu;

			}

		}

	}


	if (best_cost > cost_tabu) {

		best_cost = cost_tabu;

//		cout << " improve " << endl;
//		cout << " best cost " << best_cost << endl;

	}

	if (!S.is_feasible(I, getk())) {
		cerr << " NOT FEASIBLE!!! " << endl;
	}


	return;

}


