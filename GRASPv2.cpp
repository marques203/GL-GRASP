/*
 * GRASPv2.cpp
 *
 *  Created on: 23 mag 2017
 *      Author: antonio
 */

#include "GRASPv2.h"

GRASP_v2::GRASP_v2(HDAG &I, char *argv[]) {
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

GRASP_v2::~GRASP_v2() {
	// TODO Auto-generated destructor stub
}

void GRASP_v2::algorithm(HDAG &I, int argc, char *argv[]) {

	SolutionIMLCM s;
	set_best_cost(s.getCost(I));
	set_time_to_best(0.001);

	double initial_alpha = getalpha();

	if (argc < 12) {

		MTRand rng;
		set_best_out_file(argv[8]);
		set_complete_out_file(argv[9]);
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
	// "Execute grasp version 1";
	allocate_best_solution(I);
	auto *ES = getEliteSet();
	ES = new Elite_Set(atoi(argv[12]), atof(argv[13]));
	for (int i = 0; i < (int)getmaxIter(); i++) {
		setseed(i % 100);
		if (initial_alpha == 0) {
			setalpha(rng.randDblExc(1));
		}
		HDAG *S = new HDAG;
		run(I, *S);
		ES->add_to_Elite(*S, S->get_cost());
		S = nullptr;
	}
	double old_best = get_best_cost();
	auto *PR = getPR();
	PR = new Path_Relinking(argv[11]);
	PR->Relink(*ES, get_best_solution(), I);
	set_best_cost(get_best_solution().get_cost());
	set_end_iteration(chrono::system_clock::now());
	int elapsed = chrono::duration_cast<chrono::milliseconds>(
			get_end_iteration() - get_start_algorithm()).count();
	double current = (double) ((double) elapsed / (double) 1000);
	if (old_best > get_best_cost()) {
		set_time_to_best(current);
	}
	print_complete_data_out_file(get_best_cost(), current);
	print_best_data_out_file();
	close_out_file();

	return;

}

void GRASP_v2::run(HDAG &I, HDAG &S) {

	this->construction(I, S);

	if (strcmp(getlocalsearch(), "no") != 0) {

		this->localsearch(I, S);

	}
	else {

		SolutionIMLCM s;
		s.initialize_costsForLevels(I.getLEVELS().size());

		double best_cost = s.getCost(S);

		S.set_cost(best_cost);

		set_end_iteration(chrono::system_clock::now());

		int elapsed =
				chrono::duration_cast<chrono::milliseconds>(get_end_iteration() -
						get_start_algorithm()).count();
		double current = (double)((double)elapsed / (double)1000);

		if (best_cost < get_best_cost()) {

			set_best_cost(best_cost);
			set_time_to_best(current);

		}

		print_complete_data_out_file(best_cost, current);

		delete_all_iterative_structure(I);

	}

	return;

}

void GRASP_v2::construction(HDAG &I, HDAG &S) {

	initialize_all_iteration_structure(I);

	set_start_iteration(chrono::system_clock::now());

	MTRand rng(getseed());

//	auto &CL = getCL();
//	auto &v_maxoutdeg = getv_maxoutdeg();
//	auto &OUTDEG = getOD();

	allocate_solution(I, S);

	compute_rhomax(I, S);

//	build_maxdegree_node_structure(I, S);

	build_RCL(I, S);

	while (getTotalToAssign() > 0) {

		pair<unsigned, unsigned> p;

		p = extract_from_RCL(rng);

		add_node_to_solution(I, S, p.first, p.second, rng);

		S.increase_total_nodes(1);

		update_rho(I, S, p.first, p.second);

		update_RCL(I, S, p.first, p.second);

	}

	if (!S.is_feasible(I, getk()))
		cerr << " NOT FEASIBLE!!! " << endl;

	return;

}

void GRASP_v2::allocate_solution(HDAG &I, HDAG &S) {

	auto k = getk();
	auto &CL = getCL();
	auto &toAssign = gettoAssign();
	auto &totalToAssign = getTotalToAssign();
	auto &residual_shift = getresidual_shift();
	auto &upmost_inserted_node = getupmost_inserted_node();
//	auto &OD = getOD();

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

//				maxdeg = max(OD[l][i], maxdeg);

				toAssign[l]++;
				totalToAssign++;

			}

		}

	}

	return;
}

void GRASP_v2::add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i,
		MTRand &rng) {

	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	auto &CL = getCL();

	auto &upmost_inserted_node = getupmost_inserted_node();

	unsigned n_nodes_in_level = LEVELS_S[l].size();

	// Compute the barycenter
	unsigned bc = S.compute_barycenter(CL, I, l, i);

	if (bc > LEVELS_S[l].size() - 1) {
		bc = LEVELS_S[l].size() - 1;
	}

	bc = find_feasible_position_incremental_node(I, S, l, i, bc);

	// If the new position is free then a shifting
	// should not be performed, otherwise it is required
/*	if (!is_position_free(S, l, bc)) {

		// Find first occupied node under upmost_inserted_node[l]
		unsigned first_occupied = upmost_inserted_node[l] - 1;
		while (first_occupied >= 0 && is_position_free(S, l, first_occupied)) {
			first_occupied--;
		}

		if (first_occupied == upmost_inserted_node[l] - 1
				&& upmost_inserted_node[l] != n_nodes_in_level - 1) {
			// the first node under up-most is occupied
			// then a shift is required and the up-most does not
			// occupy the last position in the layer

			unsigned pos_upmost = upmost_inserted_node[l];
			unsigned id_upmost = IDs_S[l][pos_upmost];
			Pos_S[l][id_upmost] = pos_upmost + 1;
			LEVELS_S[l][pos_upmost + 1].swap(LEVELS_S[l][pos_upmost]);
			Os_S[l][pos_upmost + 1] = Os_S[l][pos_upmost];
			IDs_S[l][pos_upmost + 1] = id_upmost;

			// empties the swapped level
			Os_S[l][pos_upmost] = n_nodes_in_level;
			IDs_S[l][pos_upmost] = n_nodes_in_level;

			// update upmost_inserted
			upmost_inserted_node[l]++;

		} else if (Os[l][bc] == 0 && bc == upmost_inserted_node[l]
				&& upmost_inserted_node[l] != n_nodes_in_level - 1) {

			// New inserted node will be put in the up-most position
			// and the older up-most node does not occupy the last
			// position in the layer and it is not an original node

			unsigned pos_upmost = upmost_inserted_node[l];
			unsigned id_upmost = IDs_S[l][pos_upmost];
			Pos_S[l][id_upmost] = pos_upmost + 1;
			LEVELS_S[l][pos_upmost + 1].swap(LEVELS_S[l][pos_upmost]);
			Os_S[l][pos_upmost + 1] = Os_S[l][pos_upmost];
			IDs_S[l][pos_upmost + 1] = id_upmost;

			// empties the swapped level
			Os_S[l][pos_upmost] = n_nodes_in_level;
			IDs_S[l][pos_upmost] = n_nodes_in_level;

			if (pos_upmost + 1 > upmost_inserted_node[l])
				upmost_inserted_node[l] = pos_upmost + 1;

		} else if (first_occupied == upmost_inserted_node[l] - 1
				&& upmost_inserted_node[l] == n_nodes_in_level - 1) {
			unsigned last_free = first_occupied - 1;

			// up-most occupy the last position in the level
			while (last_free >= 0 && is_position_free(S, l, last_free)) {
				last_free--;
			}
			first_occupied = last_free;

		} else {

			first_occupied = upmost_inserted_node[l];
		}

		// if now bc is free can insert the new node
		// else
		if (!is_position_free(S, l, bc)) {

			for (unsigned i_ = first_occupied; i_ >= bc; --i_) {

				if (IDs_S[l][i_] != n_nodes_in_level) {

					if (Os_S[l][i_]) {

						unsigned id_i_ = IDs_S[l][i_];
						Pos_S[l][id_i_] = i_ + 1;
						LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
						Os_S[l][i_ + 1] = Os_S[l][i_];
						IDs_S[l][i_ + 1] = id_i_;

						// empties the swapped level
						Os_S[l][i_] = n_nodes_in_level;
						IDs_S[l][i_] = n_nodes_in_level;

						if (is_position_free(S, l, bc)) break;


					} else {

						unsigned finish_incremental_block = i_;
						bool space_after_incremental_block = true;

						while (finish_incremental_block > 0) {
							finish_incremental_block--;
							if (IDs_S[l][finish_incremental_block] ==
									n_nodes_in_level) {
								break;
							}
							else if (Os_S[l][finish_incremental_block] ||
									finish_incremental_block == bc) {
								space_after_incremental_block = false;
								break;
							}
						}

						if (finish_incremental_block > 0 &&
								space_after_incremental_block) {
							i_ = finish_incremental_block;
						}
						else {

							unsigned id_i_ = IDs_S[l][i_];
							Pos_S[l][id_i_] = i_ + 1;
							LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
							Os_S[l][i_ + 1] = Os_S[l][i_];
							IDs_S[l][i_ + 1] = id_i_;

							// empties the swapped level
							Os_S[l][i_] = n_nodes_in_level;
							IDs_S[l][i_] = n_nodes_in_level;

							if (is_position_free(S, l, bc)) break;

						}

					}

				}

				if (i_ == 0)
					break;

			}

		}

	}

	if (bc > upmost_inserted_node[l]) {
		upmost_inserted_node[l] = bc;
	}*/


//	cerr << " pos selected " << bc << endl;

	// If the bc position is free then a shifting
	// should not be performed, otherwise it is required
	/*if (!is_position_free(S, l, bc)) {

		// Find first occupied node under upmost_inserted_node[l]
		unsigned first_occupied = upmost_inserted_node[l] - 1;
		while (first_occupied >= 0 && is_position_free(S, l, first_occupied)) {
			first_occupied--;
		}

		if (first_occupied == upmost_inserted_node[l] - 1 &&
				upmost_inserted_node[l] != n_nodes_in_level - 1) {
			// the first node under upmost is occupied
			// then a shift is required

			unsigned pos_upmost = upmost_inserted_node[l];
			unsigned id_upmost = IDs_S[l][pos_upmost];
			Pos_S[l][id_upmost] = pos_upmost + 1;
			LEVELS_S[l][pos_upmost + 1].swap(LEVELS_S[l][pos_upmost]);
			Os_S[l][pos_upmost + 1] = Os_S[l][pos_upmost];
			IDs_S[l][pos_upmost + 1] = id_upmost;

			// empties the swapped level
			Os_S[l][pos_upmost] = n_nodes_in_level;
			IDs_S[l][pos_upmost] = n_nodes_in_level;

			// update upmost_inserted
			upmost_inserted_node[l]++;


		} else if (Os[l][bc] == 0 && bc == upmost_inserted_node[l] &&
				upmost_inserted_node[l] != n_nodes_in_level - 1) {

			unsigned pos_upmost = upmost_inserted_node[l];
			unsigned id_upmost = IDs_S[l][pos_upmost];
			Pos_S[l][id_upmost] = pos_upmost + 1;
			LEVELS_S[l][pos_upmost + 1].swap(LEVELS_S[l][pos_upmost]);
			Os_S[l][pos_upmost + 1] = Os_S[l][pos_upmost];
			IDs_S[l][pos_upmost + 1] = id_upmost;

			// empties the swapped level
			Os_S[l][pos_upmost] = n_nodes_in_level;
			IDs_S[l][pos_upmost] = n_nodes_in_level;

			if (pos_upmost + 1 > upmost_inserted_node[l])
				upmost_inserted_node[l] = pos_upmost + 1;

		} else if (first_occupied == upmost_inserted_node[l] - 1 &&
				upmost_inserted_node[l] == n_nodes_in_level - 1) {
				unsigned last_free = first_occupied-1;
				while (last_free >= 0 && !is_position_free(S, l, last_free)) {
					last_free--;
				}
				first_occupied = last_free - 1;
		}
		else {
			first_occupied = upmost_inserted_node[l];
		}

		// if now bc is free can insert the new node
		// else
		if (!is_position_free(S, l, bc)) {

			for (unsigned i_ = first_occupied; i_ >= bc; --i_) {

				if (IDs_S[l][i_] != n_nodes_in_level) {

					if (Os_S[l][i_]) {

						unsigned id_i_ = IDs_S[l][i_];
						Pos_S[l][id_i_] = i_ + 1;
						LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
						Os_S[l][i_ + 1] = Os_S[l][i_];
						IDs_S[l][i_ + 1] = id_i_;

						// empties the swapped level
						Os_S[l][i_] = n_nodes_in_level;
						IDs_S[l][i_] = n_nodes_in_level;

						//						}

					} else {

						if ((i_ > 0 && IDs_S[l][i_ - 1] != n_nodes_in_level) || i_ == bc) {

							unsigned id_i_ = IDs_S[l][i_];
							Pos_S[l][id_i_] = i_ + 1;
							LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
							Os_S[l][i_ + 1] = Os_S[l][i_];
							IDs_S[l][i_ + 1] = id_i_;

							// empties the swapped level
							Os_S[l][i_] = n_nodes_in_level;
							IDs_S[l][i_] = n_nodes_in_level;

						}

					}

				}

				if (i_ == 0) break;

			}

		}

	}*/

	if (!is_position_free(S, l, bc)) {

		// Find first free position up to new_pos
		unsigned i_ = bc;
		while (true) {
			i_++;
			if (is_position_free(S, l, i_))
				break;
		}

		unsigned first_free = i_;
		unsigned first_busy = first_free - 1;

		for (i_ = first_busy; i_ >= bc; --i_) {

			unsigned id_i_ = IDs_S[l][i_];
			Pos_S[l][id_i_] = i_ + 1;
			LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
			Os_S[l][i_ + 1] = Os_S[l][i_];
			IDs_S[l][i_ + 1] = id_i_;

			// empties the swapped level
			Os_S[l][i_] = n_nodes_in_level;
			IDs_S[l][i_] = n_nodes_in_level;

			if (is_position_free(S, l, bc)) break;

			if (i_ == 0)
				break;

		}

	}


	if (bc > upmost_inserted_node[l]) {
		upmost_inserted_node[l] = bc;
	}


	unsigned id_u = IDs[l][i];
	IDs_S[l][bc] = id_u;
	Pos_S[l][id_u] = bc;
	Os_S[l][bc] = Os[l][i];
	add_link(I, S, l, i, bc);

	// Remove from CL
	CL[l][i] = 0;

	gettoAssign()[l]--;
	getTotalToAssign()--;


//	 cerr << " node inserted! " << endl;
//	 S.printHDAG(S, "HDAG-construction-2.tex");
//	 cout << " graph printed!" << endl;
//	 getchar();

	return;

}

unsigned GRASP_v2::find_feasible_position_incremental_node(HDAG &I, HDAG &S,
		unsigned l, unsigned i, unsigned bc) {

	auto &LEVELS = I.getLEVELS();

	unsigned new_bc = bc;
	bool direction = true;
	unsigned up = 1;
	unsigned down = 1;

	if (is_position_free(S, l, bc))
		return bc;

	unsigned violation = check_violation(I, S, l, i, bc);

	while (violation != 0) {

		new_bc = bc;

		if (direction) {

			new_bc = new_bc + up;
			up++;
			direction = false;

		} else {

			new_bc = new_bc - down;
			down++;
			direction = true;

		}

		if (new_bc >= 0 && new_bc < LEVELS[l].size()) {
			violation = check_violation(I, S, l, i, new_bc);
		}

	}

	return new_bc;

}

unsigned GRASP_v2::check_violation(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned bc) {

	auto k = getk();

	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	unsigned n_nodes_in_lev = LEVELS_S[l].size();

	if (is_position_free(S, l, bc))
		return 0;

	// Otherwise, if the position is occupied
	// check if the nodes from bc to the last
	// occupied position in the level can be
	// shifted up of one position
	unsigned last_in_block = bc;
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

	/*	// Find new block
	 while (last_in_block < n_nodes_in_lev && is_position_free(S, l, last_in_block)) {
	 last_in_block++;
	 }
	 if (last_in_block != n_nodes_in_lev) {
	 first_in_block = last_in_block;
	 }
	 while (last_in_block < n_nodes_in_lev && !is_position_free(S, l, last_in_block)) {
	 last_in_block++;
	 }*/

	/*unsigned upper_bound_pos_in_level = LEVELS_S[l].size();
	 for (unsigned i_ = bc; i_ < upper_bound_pos_in_level; ++i_) {

	 unsigned id = IDs_S[l][i_];

	 if (id != upper_bound_pos_in_level) {

	 // If the node is an original node
	 // have to consider an eventual shift
	 if (Os_S[l][i_]) {

	 // The original position of the node
	 // in the the instance
	 unsigned pos = Pos[l][id];

	 // The new potential position in case
	 // of shift
	 unsigned pos_ = i_ + 1;

	 // Maximum feasible shifting
	 unsigned pos_lim = pos + k;

	 // The shifting leads to an overflow
	 if (pos_ >= upper_bound_pos_in_level) {
	 return 1;
	 }

	 // The shifting leads to a violation of
	 // the constraints position
	 if (pos_ > pos_lim) {
	 return 1;
	 }

	 }// The actual node that has to be
	 // evaluated for up shift is an
	 // incremental node
	 else {

	 // If the previous position is free
	 // then a shifting is not necessary
	 // in this case ... otherwise it must
	 // to be considered
	 if (!is_position_free(S, l, i_ - 1)) {

	 // The shifting leads to an overflow
	 if (i_ + 1 > upper_bound_pos_in_level) return 1;

	 }

	 }

	 }

	 }*/

	// The position is feasible, then the algorithm
	// can be proceeds to shifting up of one position the node
	return 0;

}

/*
 void GRASP_v2::construction(HDAG &I, HDAG &S) {

 initialize_all_iteration_structure(I);

 set_start_iteration(clock());

 MTRand rng(getseed());

 auto &CL = getCL();
 auto &v_maxoutdeg = getv_maxoutdeg();
 auto &OUTDEG = getOD();

 allocate_solution(I, S);

 compute_rhomax(I, S);

 build_maxdegree_node_structure(I, S);

 build_RCL(I, S);

 while (getTotalToAssign() > 0) {

 pair<unsigned, unsigned> p;

 if (getTotalInRCL() > 0) {

 // Select an element from the RCL
 p = extract_from_RCL(rng);

 }
 else {

 while (true) {
 if (v_maxoutdeg.size() > 0) {
 p = extract_maxdegree_node(I, rng);
 if (CL[p.first][p.second] == 1) {
 break;
 }
 }
 else {
 unsigned maximum_available_v_degree = 0;
 for (unsigned l = 0; l < OUTDEG.size(); ++l) {
 for (unsigned i = 0; i < OUTDEG[l].size(); ++i) {
 if (CL[l][i] == 1) {
 if (maximum_available_v_degree < OUTDEG[l][i]) {
 maximum_available_v_degree = OUTDEG[l][i];
 p.first = l;
 p.second = i;
 }
 }
 }
 }
 break;
 }
 }

 }

 add_node_to_solution(I, S, p.first, p.second, rng);

 update_rhomax(I, S, p.first, p.second);

 update_rho(I, S, p.first, p.second);

 }

 //	S.printHDAG(S, "HDAG-inst.tex");
 //	cout << "HEERE" << endl;
 //	getchar();

 return;

 }


 void GRASP_v2::allocate_solution(HDAG &I, HDAG &S) {

 auto k = getk();
 auto &CL = getCL();
 auto &maxoutdeg = getmaxoutdeg();
 auto &toAssing = gettoAssign();
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

 if (upmost_inserted_node[l] == vn) {
 upmost_inserted_node[l] = i;
 }
 else if (upmost_inserted_node[l] < i) {
 upmost_inserted_node[l] = i;
 }

 }
 else {

 maxoutdeg = max(OD[l][i], maxoutdeg);

 toAssing[l]++;
 totalToAssign++;

 }

 }

 }

 }


 void GRASP_v2::add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng) {

 auto &LEVELS = I.getLEVELS();
 auto &IDs = I.getIDs();
 auto &Os = I.getOs();

 auto &IDs_S = S.getIDs();
 auto &Pos_S = S.getPos();
 auto &Os_S = S.getOs();

 auto &CL = getCL();

 auto &upmost_inserted_node = getupmost_inserted_node();

 // Compute the barycenter
 unsigned bc = S.compute_barycenter(l);

 bc = find_feasible_position_incremental_node(I, S, l, i, bc);

 cout << " total to assign " << getTotalToAssign() << endl;

 unsigned old_upmost = upmost_inserted_node[l];

 // If the bc position is free then a shifting
 // should not be performed, otherwise it is required
 if (!is_position_free(S, l, bc)) {

 S.printHDAG(S, "HDAG-inst.tex");
 cout << " i will insert the node " << I.getIDs()[l][i] << " in the lev " << l << endl;
 cout << " it position will be " << bc << endl;
 cout << " upmost_inserted_node in " << l << " " << upmost_inserted_node[l] << endl;
 getchar();
 for (unsigned i_ = upmost_inserted_node[l]; i_ >= bc; i_--) {

 if (IDs_S[l][i_] != LEVELS[l].size() && !S.isOriginalNode(i_, l)) {
 if (!is_position_free(S, l, i_ - 1)) {
 S.swapPositions(i_, i_+1, l);
 if (I.getIDs()[l][i] == 35) {
 cout << " swapping incremental concrete " << endl;
 S.printHDAG(S, "HDAG-inst.tex");
 getchar();
 }
 }
 else{
 if (I.getIDs()[l][i] == 35) {
 cout << " swapping incremental ficticiuos " << endl;
 S.printHDAG(S, "HDAG-inst.tex");
 getchar();
 }
 }
 }
 else {
 S.swapPositions(i_, i_+1, l);
 if (I.getIDs()[l][i] == 35) {
 cout << " swapping original " << endl;
 S.printHDAG(S, "HDAG-inst.tex");
 getchar();
 }
 }
 }
 }

 if (upmost_inserted_node[l] == LEVELS[l].size()) {
 upmost_inserted_node[l] = bc;
 }
 else if (upmost_inserted_node[l] < bc) {
 upmost_inserted_node[l] = bc;
 }
 else if (upmost_inserted_node[l] < LEVELS[l].size() - 1
 && IDs_S[l][upmost_inserted_node[l]+1] != LEVELS[l].size()) {
 upmost_inserted_node[l]++;
 }

 unsigned id_u = IDs[l][i];
 IDs_S[l][bc] = id_u;
 Pos_S[l][id_u] = bc;
 Os_S[l][bc] = Os[l][i];
 add_link(I, S, l, i, bc);

 // Remove from CL
 CL[l][i] = 0;

 gettoAssign()[l]--;
 getTotalToAssign()--;

 return;

 }

 unsigned GRASP_v2::is_feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) {

 auto k = getk();

 auto &Pos = I.getPos();

 auto &LEVELS_S = S.getLEVELS();
 auto &IDs_S = S.getIDs();
 auto &Os_S = S.getOs();

 if (is_position_free(S, l, bc)) return 0;

 // Otherwise, if the position is occupied
 // check if the nodes from bc to the last
 // occupied position in the level can be
 // shifted up of one position
 unsigned upper_bound_pos_in_level = LEVELS_S[l].size();
 for (unsigned i_ = bc; i_ < upper_bound_pos_in_level; ++i_) {

 unsigned id = IDs_S[l][i_];

 if (id != upper_bound_pos_in_level) {

 // If the node is an original node
 // have to consider an eventual shift
 if (Os_S[l][i_]) {

 // The original position of the node
 // in the the instance
 unsigned pos = Pos[l][id];

 // The new potential position in case
 // of shift
 unsigned pos_ = i_ + 1;

 // Maximum feasible shifting
 unsigned pos_lim = pos + k;

 // The shifting leads to an overflow
 if (pos_ >= upper_bound_pos_in_level) {
 return 1;
 }

 // The shifting leads to a violation of
 // the constraints position
 if (pos_ > pos_lim) {
 return 1;
 }

 }// The actual node that has to be
 // evaluated for up shift is an
 // incremental node
 else {

 // If the previous position is free
 // then a shifting is not necessary
 // in this case ... otherwise it must
 // to be considered
 if (!is_position_free(S, l, i_ - 1)) {

 // The shifting leads to an overflow
 if (i_ + 1 > upper_bound_pos_in_level) return 1;

 }
 }

 }

 }

 // The position is feasible, then the algorithm
 // can be proceeds to shifting up of one position the node
 return 0;


 }

 unsigned GRASP_v2::find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) {

 auto &LEVELS = I.getLEVELS();

 unsigned new_bc = bc;
 bool direction = true;
 unsigned up = 1;
 unsigned down = 1;

 if (is_position_free(S, l, bc))
 return bc;

 unsigned violation = is_feasible_position(I, S, l, i, bc);

 while (violation != 0) {

 new_bc = bc;

 if (direction) {

 new_bc = new_bc + up;
 up++;
 direction = false;

 } else {

 new_bc = new_bc - down;
 down++;
 direction = true;

 }

 if (new_bc >= 0 && new_bc < LEVELS[l].size()) {
 violation = is_feasible_position(I, S, l, i, new_bc);
 }

 }

 return new_bc;


 }
 */
