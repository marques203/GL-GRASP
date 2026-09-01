/*
 * GRASPv3.cpp
 *
 *  Created on: 11 lug 2017
 *      Author: antonio
 */

#include "GRASPv3.h"

GRASP_v3::GRASP_v3(HDAG &I, char *argv[]) :
tau(0), global_min_g(INT32_MAX), global_max_g(0) {
	// TODO Auto-generated constructor stub

	setk(atoi(argv[3]));
	setseed(1);
	setalpha(atof(argv[4]));
	setmaxIter(atoi(argv[6]));
	setlocalsearch(argv[5]);
	setTotalInRCL(0);
	setTotalToAssign(0);
	set_best_cost(INT64_MAX);
	set_start_algorithm(chrono::system_clock::now());
	allocate_and_build_all_structure(I);
	compute_degree_matrix(I);

}

void GRASP_v3::algorithm(HDAG &I, int argc, char *argv[]) {

	SolutionIMLCM s;
	set_best_cost(s.getCost(I));
	set_time_to_best(0.001);
	cout << "instance cost " << get_best_cost();

	double initial_alpha = getalpha();

	if (argc < 12) {

		MTRand rng;
		set_best_out_file(argv[8]);
		set_complete_out_file(argv[9]);
		setTimeLimit(atof(argv[10]));
		// "Execute grasp version 1";
		allocate_best_solution(I);
		for (int i = 0; i < (int)getmaxIter(); i++) {
			setseed(i % 100);
			if (initial_alpha == 0) {
				setalpha(rng.randDblExc(1));
			}
			HDAG *S = new HDAG;
			run(I, *S);
			S = nullptr;
			int elapsed = chrono::duration_cast<chrono::milliseconds>(
					get_end_iteration() - get_start_algorithm()).count();
			double current = (double) ((double) elapsed / (double) 1000);
			if (current >= getTimeLimit()) {
				break;
			}
		}

		int elapsed = chrono::duration_cast<chrono::milliseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000);
		print_complete_data_out_file(get_best_cost(), current);
		print_best_data_out_file();
		close_out_file();

		return;

	}

	MTRand rng;
	set_best_out_file(argv[8]);
	set_complete_out_file(argv[9]);
	setTimeLimit(atof(argv[10]));
	// "Execute grasp version 1";
	allocate_best_solution(I);
	auto *ES = getEliteSet();
	ES = new Elite_Set(atoi(argv[12]), atof(argv[13]));

//	chrono::time_point<chrono::system_clock> init_cLS = chrono::system_clock::now();

	for (int i = 0; i < (int)getmaxIter(); i++) {
		setseed(i % 100);
		if (initial_alpha == 0) {
			setalpha(rng.randDblExc(1));
		}
		HDAG *S = new HDAG;
		run(I, *S);
		ES->add_to_Elite(*S, S->get_cost());
		S = nullptr;
		double elapsed = chrono::duration_cast<chrono::milliseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000);
		if (current > getTimeLimit()) {
			break;
		}
	}

//	chrono::time_point<chrono::system_clock> end_cLS = chrono::system_clock::now();
//	chrono::time_point<chrono::system_clock> init_PR = chrono::system_clock::now();

	double old_best = get_best_cost();
	auto *PR = getPR();
	PR = new Path_Relinking(argv[11]);
	PR->Relink(*ES, get_best_solution(), I);
	set_best_cost(get_best_solution().get_cost());
	set_end_iteration(chrono::system_clock::now());
	double elapsed = chrono::duration_cast<chrono::milliseconds>(
			get_end_iteration() - get_start_algorithm()).count();
	double current = (double) ((double) elapsed / (double) 1000);
	if (old_best > get_best_cost()) {
		set_time_to_best(current);
	}
	print_complete_data_out_file(get_best_cost(), current);
	print_best_data_out_file();
	close_out_file();

/*	chrono::time_point<chrono::system_clock> end_PR = chrono::system_clock::now();

	double totDur = chrono::duration_cast<chrono::milliseconds>(
			end_PR - init_cLS).count();
	totDur = totDur / 1000;

	double cLSDur = chrono::duration_cast<chrono::milliseconds>(
			end_cLS - init_cLS).count();
	cLSDur = cLSDur / 1000;

	double PRDur = chrono::duration_cast<chrono::milliseconds>(
			end_PR - init_PR).count();
	PRDur = PRDur / 1000;

	FILE *F = fopen("GRASP3PRtimeDistribution.csv", "a");

	fprintf(F, "%lf, %lf, %lf\n", cLSDur, PRDur, totDur);

	fclose(F);*/

	return;

}

void GRASP_v3::run(HDAG &I, HDAG &S) {


	this->construction(I, S);

	if (strcmp(getlocalsearch(), "no") != 0) {

		this->localsearch(I, S);


	} else {

		SolutionIMLCM s;
		s.initialize_costsForLevels(I.getLEVELS().size());

		double best_cost = s.getCost(S);

		S.set_cost(best_cost);

		set_end_iteration(chrono::system_clock::now());

		int elapsed = chrono::duration_cast<chrono::milliseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000);

		if (best_cost < get_best_cost()) {

			set_best_cost(best_cost);
			set_time_to_best(current);

		}


		print_complete_data_out_file(best_cost, current);

		delete_all_iterative_structure(I);

	}

	return;

}

void GRASP_v3::construction(HDAG &I, HDAG &S) {

	initialize_all_iteration_structure(I);

	set_start_iteration(chrono::system_clock::now());

	MTRand rng(getseed());

	allocate_solution(I, S);

	build_UNASSIGNED(I, S);

	updateTau();

	build_RCL_for_GRASPv3(I, S);

	while (getTotalToAssign() > 0) {

		pair<unsigned, unsigned> p;
		p = extract_from_RCL(rng);

		unsigned new_pos = find_new_position(p.first, p.second);

		add_node_to_solution_for_GRASPv3(I, S, p.first, p.second, new_pos);

		apply_updates_after_insertion(I, S, p.first, p.second);

		double old_tau = tau;
		updateTau();

		re_build_RCL_for_GRASPv3(I, S, old_tau, p.first, p.second);

	}

	free_UNASSIGNED();

	if (!S.is_feasible(I, getk()))
		cerr << " NOT FEASIBLE!!! " << endl;

	return;

}

void GRASP_v3::allocate_solution(HDAG &I, HDAG &S) {

	auto k = getk();
	auto &CL = getCL();
	auto &maxdeg = getmaxdeg();
	auto &toAssign = gettoAssign();
	auto &totalToAssign = getTotalToAssign();
	auto &residual_shift = getresidual_shift();
	auto &upmost_inserted_node = getupmost_inserted_node();
	auto &OD = getOD();

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	unsigned ln, vn;

	ln = I.getLevNumber();

	S.allocateLEVELS(ln);
	S.allocateIDs(ln);
	S.allocateOs(ln);
	S.allocatePos(ln);

	for (unsigned l = 0; l < ln; ++l) {

		vn = LEVELS[l].size();

		S.allocateLevel(l, vn);
		S.allocateIDs(l, vn);
		S.allocateOs(l, vn);
		S.allocatePos(l, vn);
		residual_shift[l] = k;
		upmost_inserted_node[l] = vn;

		for (unsigned i = 0; i < vn; ++i) {

			if (I.isOriginalNode(i, l)) {

				unsigned id_u = IDs[l][i];
				IDs_S[l][i] = id_u;
				Pos_S[l][id_u] = i;
				Os_S[l][i] = Os[l][i];
				add_link(I, S, l, i, i);

				// Remove from CL
				CL[l][i] = 0;

				S.increase_total_nodes(1);

				if (upmost_inserted_node[l] == vn) {
					upmost_inserted_node[l] = i;
				} else if (upmost_inserted_node[l] < i) {
					upmost_inserted_node[l] = i;
				}

			} else {

				maxdeg = max(OD[l][i], maxdeg);

				toAssign[l]++;
				totalToAssign++;

			}

		}

	}

	return;

}

void GRASP_v3::add_node_to_solution_for_GRASPv3(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned new_pos) {

	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	auto &CL = getCL();

	auto &upmost_inserted_node = getupmost_inserted_node();

	unsigned n_nodes_in_level = LEVELS_S[l].size();

	if (!is_position_free(S, l, new_pos)) {

		// Find first free position up to new_pos
		unsigned i_ = new_pos;
		while (true) {
			i_++;
			if (is_position_free(S, l, i_))
				break;
		}

		unsigned first_free = i_;
		unsigned first_busy = first_free - 1;

		for (i_ = first_busy; i_ >= new_pos; --i_) {

			unsigned id_i_ = IDs_S[l][i_];
			Pos_S[l][id_i_] = i_ + 1;
			LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
			Os_S[l][i_ + 1] = Os_S[l][i_];
			IDs_S[l][i_ + 1] = id_i_;

			// empties the swapped level
			Os_S[l][i_] = n_nodes_in_level;
			IDs_S[l][i_] = n_nodes_in_level;

			if (is_position_free(S, l, new_pos)) break;

			if (i_ == 0)
				break;

		}

	}

	if (new_pos > upmost_inserted_node[l]) {
		upmost_inserted_node[l] = new_pos;
	}

	unsigned id_u = IDs[l][i];
	IDs_S[l][new_pos] = id_u;
	Pos_S[l][id_u] = new_pos;
	Os_S[l][new_pos] = Os[l][i];
	add_link(I, S, l, i, new_pos);

	// Remove from CL
	CL[l][i] = 0;

	gettoAssign()[l]--;
	getTotalToAssign()--;


	return;

}

unsigned GRASP_v3::check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i,
		unsigned pos) {

	auto k = getk();

	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	unsigned n_nodes_in_lev = LEVELS_S[l].size();

	if (is_position_free(S, l, pos))
		return 0;

	// check if the nodes from bc to the last
	// occupied position in the level can be
	// shifted up of one position
	unsigned last_in_block = pos;
	while (last_in_block < n_nodes_in_lev
			&& !is_position_free(S, l, last_in_block)) {
		if (Os_S[l][last_in_block]) {

			unsigned id_ = IDs_S[l][last_in_block];
			unsigned origin_pos = Pos[l][id_];
			unsigned potential_new_pos = last_in_block + 1;
			unsigned max_feasibile_pos = origin_pos + k;

			if (potential_new_pos > max_feasibile_pos)
				return 1;

		}
		last_in_block++;
	}

	// The whole block cannot be shifted up
	if (last_in_block == n_nodes_in_lev)
		return 1;

	// The position is feasible, then the algorithm
	// can be proceeds to shifting up of one position the node
	return 0;

}

void GRASP_v3::build_UNASSIGNED(HDAG &I, HDAG &S) {

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &CL = getCL();

	auto &LEVELS_S = S.getLEVELS();

	UNASSIGNED.resize(LEVELS.size());

	min_g.resize(LEVELS.size());

	for (unsigned l = 0; l < LEVELS.size(); ++l) {

		UNASSIGNED[l].resize(LEVELS[l].size());
		min_g[l].resize(LEVELS[l].size(), INT32_MAX);

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			unsigned id = IDs[l][i];

			unsigned original_p = Pos[l][id];

			if (CL[l][original_p] == 1) {

				// Consider all the potential positions for the node id
				for (unsigned i_ = 0; i_ < LEVELS_S[l].size(); ++i_) {

					unsigned violation = check_violation(I, S, l, original_p,
							i_);

					if (!violation) {

						double potential_cost =
								compute_Bward_potential_cost(I, S, l, original_p, i_)
							  + compute_Fward_potential_cost(I, S, l, original_p, i_);

						tuple<unsigned, unsigned, double, bool> t;
						get<0>(t) = id;
						get<1>(t) = i_;
						get<2>(t) = potential_cost;
						get<3>(t) = 1;

						if (potential_cost < min_g[l][i]) {
							min_g[l][i] = potential_cost;
						}

						UNASSIGNED[l][i].push_back(t);

						if (potential_cost < global_min_g) {
							global_min_g = potential_cost;
							//position_min_g = make_tuple(l, i,
							//		UNASSIGNED[l][i].size() - 1);
							get<0>(position_min_g) = l;
							get<1>(position_min_g) = i;
							get<2>(position_min_g) = UNASSIGNED[l][i].size() - 1;
						}

						if (min_g[l][i] > global_max_g) {
							global_max_g = min_g[l][i];
							//position_max_g = make_tuple(l, i,
							//		UNASSIGNED[l][i].size() - 1);
							get<0>(position_max_g) = l;
							get<1>(position_max_g) = i;
							get<2>(position_max_g) = UNASSIGNED[l][i].size() - 1;

						}

					}

				}

			}

		}

	}

	return;

}

double GRASP_v3::compute_Bward_potential_cost(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned new_p) {

	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();

	if (l == 0)
		return 0;

	unsigned nodes_in_previous_level = LEVELS_S[l - 1].size();

	auto &CL = getCL();

	double cost = 0;

	// indicate the node in the original position (l , i)
	// with "j1"

	for (unsigned j = 0; j < B_LEVELS[l][i].size(); ++j) {
		// indicate the node in the j-th position of
		// the backward star of j1 with "i1"
		unsigned i1 = B_LEVELS[l][i][j];
		unsigned original_pos_i1 = Pos[l - 1][i1];

		// check if the node i1 is in the solution
		if (CL[l - 1][original_pos_i1] == 0) {

			// now consider the position of i1
			// in the new solution "S"
			unsigned new_p_i1 = Pos_S[l - 1][i1];

			// the "potential" position of j1
			// in the new solution is new_p
			unsigned new_p_j1 = new_p;

			// now analyze all the edges (i2, j2)
			// in the new solution e checks if
			// the edges (i1, j1) and (i2, j2)
			// are "potentially" crossings

			for (unsigned i_ = 0; i_ < LEVELS_S[l - 1].size(); ++i_) {

				// indicate with i2 the node in position i_
				// in the new solution
				unsigned i2 = IDs_S[l - 1][i_];

				// check if the node is already in solution
				if (i2 != nodes_in_previous_level) {

					// consider its forward star
					// if i2 is not the same node
					// indicated by i1
					if (i2 != i1) {

						// indicate with new_p_i2 the position
						// of the node i2 in the new solution
						unsigned new_p_i2 = Pos_S[l - 1][i2];

						for (unsigned j_ = 0; j_ < LEVELS_S[l - 1][i_].size();
								++j_) {

							// indicate with j2 the j_-th node in the forward list
							// of i2
							unsigned j2 = LEVELS_S[l - 1][i_][j_];

							// indicate with new_p_j2 the position
							// of the node j2 in the new solution
							unsigned new_p_j2 = Pos_S[l][j2];

							// now check if the edges (i1, j1) and (i2, j2)
							// are crossings edges
							if (((new_p_i1 < new_p_i2) && (new_p_j1 > new_p_j2))
									|| ((new_p_i1 > new_p_i2)
											&& (new_p_j1 < new_p_j2))) {
								cost = cost + 1;

							}

						}

					}

				}

			}

		}

	}

	return cost;

}

double GRASP_v3::compute_Fward_potential_cost(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned new_p) {

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();

	unsigned nodes_in_actual_level = LEVELS_S[l].size();

	auto &CL = getCL();

	double cost = 0;

	// indicate the node in the original position (l , i)
	// with "i1"
	unsigned i1 = IDs[l][i];

	// Analyze its forward star
	for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

		// Indicate with j1 the j-th node in the
		// forward star of i1
		unsigned j1 = LEVELS[l][i][j];

		// Find its original position
		unsigned original_p_j1 = Pos[l + 1][j1];

		// Check if j1 is in the new solution
		if (CL[l + 1][original_p_j1] == 0) {

			// now consider the position of i1
			// in the new solution "S"
			unsigned new_p_i1 = new_p;

			// the "potential" position of j1
			// in the new solution "S"
			unsigned new_p_j1 = Pos_S[l + 1][j1];

			// now analyze all the edges (i2, j2)
			// in the new solution e checks if
			// the edges (i1, j1) and (i2, j2)
			// are "potentially" crossings

			for (unsigned i_ = 0; i_ < LEVELS_S[l].size(); ++i_) {

				// indicate with i2 the node in position i_
				// in the new solution
				unsigned i2 = IDs_S[l][i_];

				// check if the node is already in solution
				if (i2 != nodes_in_actual_level) {

					// consider its forward star
					// if i2 is not the same node
					// indicated by i1
					if (i2 != i1) {

						// indicate with new_p_i2 the position
						// of the node i2 in the new solution
						unsigned new_p_i2 = Pos_S[l][i2];

						for (unsigned j_ = 0; j_ < LEVELS_S[l][i_].size();
								++j_) {

							// indicate with j2 the j_-th node in the forward list
							// of i2
							unsigned j2 = LEVELS_S[l][i_][j_];

							// indicate with new_p_j2 the position
							// of the node j2 in the new solution
							unsigned new_p_j2 = Pos_S[l + 1][j2];

							// now check if the edges (i1, j1) and (i2, j2)
							// are crossings edges
							if (((new_p_i1 < new_p_i2) && (new_p_j1 > new_p_j2))
									|| ((new_p_i1 > new_p_i2)
											&& (new_p_j1 < new_p_j2))) {
								cost = cost + 1;

							}

						}

					}

				}

			}

		}

	}

	return cost;

}

void GRASP_v3::build_RCL_for_GRASPv3(HDAG &I, HDAG &S) {

	auto &Pos = I.getPos();

	for (unsigned l = 0; l < UNASSIGNED.size(); ++l) {
		for (unsigned i = 0; i < UNASSIGNED[l].size(); ++i) {
			if (min_g[l][i] <= (floor(tau))) {
				unsigned id = get<0>(UNASSIGNED[l][i][0]);
				unsigned original_pos = Pos[l][id];
				add_in_RCL(l, original_pos);
			}
		}
	}

	return;

}

unsigned GRASP_v3::find_new_position(unsigned l, unsigned i) {

	unsigned new_pos=0;

	for (unsigned j = 0; j < UNASSIGNED[l][i].size(); ++j) {
		if (get<2>(UNASSIGNED[l][i][j]) == min_g[l][i] &&
			get<3>(UNASSIGNED[l][i][j]) == 1) {
			new_pos = get<1>(UNASSIGNED[l][i][j]);
			break;
		}
	}

	return new_pos;

}

void GRASP_v3::apply_updates_after_insertion(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &IDs = I.getIDs();

	auto &LEVELS_S = S.getLEVELS();

	auto &CL = getCL();

	// THE INFEASIBILITY MUST NOT BE CHANGED!
	if (l > 0) {

		for (unsigned u = 0; u < LEVELS_S[l - 1].size(); ++u) {

			if (CL[l - 1][u] == 1) {

				unsigned id = IDs[l - 1][u];

				min_g[l - 1][u] = INT32_MAX;

				for (unsigned j = 0; j < UNASSIGNED[l - 1][u].size(); ++j) {

					unsigned is_feasibile = get<3>(UNASSIGNED[l - 1][u][j]);

					if (is_feasibile) {

						unsigned pos_to_check = get<1>(UNASSIGNED[l - 1][u][j]);

						double potential_cost =
								compute_Bward_potential_cost(I, S, l - 1, u, pos_to_check)
								+ compute_Fward_potential_cost(I, S, l - 1, u, pos_to_check);

						if (potential_cost < min_g[l - 1][u]) {
							min_g[l - 1][u] = potential_cost;
						}

						get<0>(UNASSIGNED[l - 1][u][j]) = id;
						get<1>(UNASSIGNED[l - 1][u][j]) = pos_to_check;
						get<2>(UNASSIGNED[l - 1][u][j]) = potential_cost;
						get<3>(UNASSIGNED[l - 1][u][j]) = 1;

					}
				}

			}
		}
	}
	// THE INFEASIBILITY MUST NOT BE CHANGED!
	if (l < LEVELS_S.size() - 1) {

		for (unsigned u = 0; u < LEVELS_S[l + 1].size(); ++u) {

			if (CL[l + 1][u] == 1) {

				unsigned id = IDs[l + 1][u];

				min_g[l + 1][u] = INT32_MAX;

				for (unsigned j = 0; j < UNASSIGNED[l + 1][u].size(); ++j) {

					unsigned is_feasibile = get<3>(UNASSIGNED[l + 1][u][j]);

					if (is_feasibile) {

						unsigned pos_to_check = get<1>(UNASSIGNED[l + 1][u][j]);

						double potential_cost =
								compute_Bward_potential_cost(I, S, l + 1, u, pos_to_check)
								+ compute_Fward_potential_cost(I, S, l + 1, u, pos_to_check);

						if (potential_cost < min_g[l + 1][u]) {
							min_g[l + 1][u] = potential_cost;
						}

						get<0>(UNASSIGNED[l + 1][u][j]) = id;
						get<1>(UNASSIGNED[l + 1][u][j]) = pos_to_check;
						get<2>(UNASSIGNED[l + 1][u][j]) = potential_cost;
						get<3>(UNASSIGNED[l + 1][u][j]) = 1;

					}
				}

			}
		}
	}


	// THE INFEASIBILITY SHOULD BE CHANGED!
	for (unsigned u = 0; u < LEVELS_S[l].size(); ++u) {

		if (CL[l][u] == 1) {

			unsigned id = IDs[l][u];

			min_g[l][u] = INT32_MAX;

			for (unsigned j = 0; j < UNASSIGNED[l][u].size(); ++j) {

				unsigned is_feasibile = get<3>(UNASSIGNED[l][u][j]);

				if (is_feasibile) {

					unsigned pos_to_check = get<1>(UNASSIGNED[l][u][j]);

					unsigned violation = check_violation(I, S, l, u, pos_to_check);

					if (!violation) {

						double potential_cost =
							compute_Bward_potential_cost(I, S, l, u, pos_to_check)
							+ compute_Fward_potential_cost(I, S, l, u, pos_to_check);


						if (potential_cost < min_g[l][u]) {
							min_g[l][u] = potential_cost;
						}

						get<0>(UNASSIGNED[l][u][j]) = id;
						get<1>(UNASSIGNED[l][u][j]) = pos_to_check;
						get<2>(UNASSIGNED[l][u][j]) = potential_cost;
						get<3>(UNASSIGNED[l][u][j]) = 1;

					}
					else {

						get<0>(UNASSIGNED[l][u][j]) = id;
						get<1>(UNASSIGNED[l][u][j]) = pos_to_check;
						get<2>(UNASSIGNED[l][u][j]) = INT32_MAX;
						get<3>(UNASSIGNED[l][u][j]) = 0;

					}
				}
			}
		}
	}

	global_min_g = INT32_MAX;
	global_max_g = 0;

		for (unsigned l_ = 0; l_ < UNASSIGNED.size(); ++l_) {
			for (unsigned i_ = 0; i_ < UNASSIGNED[l_].size(); ++i_) {
				if (CL[l_][i_] == 1) {
					for (unsigned j = 0; j < UNASSIGNED[l_][i_].size(); ++j) {
						bool is_feasibile = get<3>(UNASSIGNED[l_][i_][j]);
						if (is_feasibile) {
							unsigned potential_cost = get<2>(UNASSIGNED[l_][i_][j]);

							if (potential_cost < global_min_g) {
								global_min_g = potential_cost;
								//position_min_g = make_tuple(l_, i_,
								//		UNASSIGNED[l_][i_].size() - 1);
								get<0>(position_min_g) = l_;
								get<1>(position_min_g) = i_;
								get<2>(position_min_g) = UNASSIGNED[l_][i_].size() - 1;
							}

							if (min_g[l_][i_] > global_max_g) {
								global_max_g = potential_cost;
								//position_max_g = make_tuple(l_, i_,
								//		UNASSIGNED[l_][i_].size() - 1);
								get<0>(position_max_g) = l_;
								get<1>(position_max_g) = i_;
								get<2>(position_max_g) = UNASSIGNED[l_][i_].size() - 1;
							}
						}
					}
				}
			}
		}

	return;

}

void GRASP_v3::re_build_RCL_for_GRASPv3(HDAG &I, HDAG &S, double  old_tau, unsigned l, unsigned i) {

	auto &CL = getCL();

	if (tau > old_tau) {

		// In this case MUST be considered ALL the nodes
		for (unsigned l_ = 0; l_ < UNASSIGNED.size(); ++l_) {
			for (unsigned i_ = 0; i_ < UNASSIGNED[l_].size(); ++i_) {

				if (CL[l_][i_] == 1) {

					for (unsigned j = 0; j < UNASSIGNED[l_][i_].size(); ++j) {
						bool feasibile = get<3>(UNASSIGNED[l_][i_][j]);

						if (feasibile) {
							if (!is_already_in_RCL(l_, i_)) {
								if (min_g[l_][i_] <= floor(tau)) {
									add_in_RCL(l_, i_);
								}
							}
						}
					}
				}
			}
		}
	}
	else {
		// Consider only the levels interested by the changes
		if (l > 0) {
			for (unsigned i_ = 0; i_ < UNASSIGNED[l - 1].size(); ++i_) {
				if (CL[l - 1][i_] == 1) {
					for (unsigned j = 0; j < UNASSIGNED[l - 1][i_].size(); ++j) {
						bool feasibile = get<3>(UNASSIGNED[l - 1][i_][j]);
						if (feasibile) {
							if (!is_already_in_RCL(l - 1, i_)) {
								if (min_g[l - 1][i_] <= floor(tau)) {
									add_in_RCL(l - 1, i_);
								}
							}
						}
					}
				}
			}
		}
		for (unsigned i_ = 0; i_ < UNASSIGNED[l].size(); ++i_) {
			if (CL[l][i_] == 1) {
				for (unsigned j = 0; j < UNASSIGNED[l][i_].size(); ++j) {
					bool feasibile = get<3>(UNASSIGNED[l][i_][j]);
					if (feasibile) {
						if (!is_already_in_RCL(l, i_)) {
							if (min_g[l][i_] <= floor(tau)) {
								add_in_RCL(l, i_);
							}
						}
					}
				}
			}
		}
	}

	return;

}

GRASP_v3::~GRASP_v3() {
	// TODO Auto-generated destructor stub
}

