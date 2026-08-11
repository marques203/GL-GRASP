/*
 * GRASPv1.cpp
 *
 *  Created on: 24 mag 2017
 *      Author: antonio
 */

#include "GRASPv1.h"

GRASP_v1::GRASP_v1(HDAG &I, char *argv[]) {
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

GRASP_v1::~GRASP_v1() {
	// TODO Auto-generated destructor stub
}

void GRASP_v1::algorithm(HDAG &I, int argc, char *argv[]) {

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

void GRASP_v1::run(HDAG &I, HDAG &S) {

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


void GRASP_v1::construction(HDAG &I, HDAG &S) {

	initialize_all_iteration_structure(I);

	set_start_iteration(chrono::system_clock::now());

	allocate_solution(I, S);

	build_maxdegree_node_structure(I, S);

	MTRand rng(getseed());

	pair<unsigned, unsigned> v = extract_maxdegree_node(I, rng);

	add_first_node_to_solution(I, S, v.first, v.second, rng);

	S.increase_total_nodes(1);

	update_rho(I, S, v.first, v.second);

//	update_rhomax(I, S, v.first, v.second);

	while (getTotalToAssign() > 0) {

//	   cerr << " RCL-size = " << getTotalInRCL() << endl;

	   pair<unsigned, unsigned> p;

		p = extract_from_RCL(rng);

		add_node_to_solution(I, S, p.first, p.second, rng);

		S.increase_total_nodes(1);

//		update_rhomax(I, S, p.first, p.second);

		update_rho(I, S, p.first, p.second);

		update_RCL(I, S, p.first, p.second);

	}

	if (!S.is_feasible(I, getk()))
		cerr << " NOT FEASIBLE!!! " << endl;

	return;

}

void GRASP_v1::allocate_solution(HDAG &I, HDAG &S) {

	auto &CL = getCL();

	auto &toAssing = gettoAssign();
	auto &totalToAssign = getTotalToAssign();

	const auto &LEVELS = I.getLEVELS();
	auto &OD = getOD();

	unsigned ln, vn, deg;

	ln = I.getLevNumber();
	upmost_original.resize(ln);

	original_to_assign.resize(ln);

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
		upmost_original[l] = vn;
		original_to_assign[l] = 0;

		for (unsigned i = 0; i < vn; ++i) {

			CL[l][i] = 1;
			toAssing[l]++;
			totalToAssign++;

			deg = I.compute_degree(l, i);

			OD[l][i] = deg;

			if (deg == I.getMaxDegree()) {

				pair<unsigned, unsigned> tmp;
				tmp.first = l;
				tmp.second = i;

			}

			if (I.isOriginalNode(i, l)) {
				if (upmost_original[l] == vn) {
					upmost_original[l] = i;
				}
				else if (upmost_original[l] < i) {
					upmost_original[l] = i;
				}
			} else
				original_to_assign[l]++;

		}

	}

	return;

}

void GRASP_v1::add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i,
		MTRand &rng) {

	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	auto &CL = getCL();

	// Compute the barycenter
	unsigned bc = S.compute_barycenter(CL, I, l, i);

	if (bc > LEVELS[l].size() - 1) {
		bc = LEVELS[l].size() - 1;
	}

	bc = find_feasible_position_node(I, S, l, i, bc);

	unsigned id_u = IDs[l][i];
	IDs_S[l][bc] = id_u;
	Pos_S[l][id_u] = bc;
	Os_S[l][bc] = Os[l][i];
	if (I.isOriginalNode(i, l)) {

		original_to_assign[l]--;

	}
	else {

	}
	add_link(I, S, l, i, bc);

	// Remove from CL
	CL[l][i] = 0;

	gettoAssign()[l]--;
	getTotalToAssign()--;

	return;

}

unsigned GRASP_v1::check_violation(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned bc) {

	auto k = getk();
	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	unsigned id_u = IDs[l][i];

	// The node is an original node
	if (I.isOriginalNode(i, l)) {
		// bc < p(u)
		if (bc < Pos[l][id_u])
			return 1;
		// bc > p(u) + k
		if (bc > Pos[l][id_u] + k)
			return 2;
		// distance from p'(u) and last position does not consent to insert
		// all the remaining original nodes v such that p(v) > p(u)
		unsigned difference_by_upmost = upmost_original[l] - Pos[l][IDs[l][i]];
		int upper_limit = (LEVELS[l].size() - 1) - difference_by_upmost;
		if ((int) bc > upper_limit)
			return 3;
		// Exists v | p'(v) < bc and p(v) > p(u)
		if (bc > 0) {
			for (unsigned i_ = bc - 1; i_ > 0; --i_) {
				if (IDs_S[l][i_] != LEVELS_S[l].size()) {
					if (Os[l][IDs_S[l][i_]]) {
						if (Pos[l][IDs_S[l][i_]] > Pos[l][id_u]) {
							return 4;
						}
					}
				}
			}
		}
		// Exists v | p'(v) > bc and p(v) < p(u)
		for (unsigned i_ = bc + 1; i_ < LEVELS_S[l].size(); ++i_) {
			if (IDs_S[l][i_] != LEVELS_S[l].size()) {
				if (Os[l][IDs_S[l][i_]]) {
					if (Pos[l][IDs_S[l][i_]] < Pos[l][id_u]) {
						return 5;
					}
				}
			}
		}
		// Exists v | p'(v) < bc and p'(v) - p(v) > p'(u) - p(u)
		int diff_u = bc - Pos[l][id_u];
		if (bc > 0) {
			for (unsigned i_ = bc - 1; i_ > 0; --i_) {
				if (IDs_S[l][i_] != LEVELS_S[l].size()) {
					unsigned id_v = IDs_S[l][i_];
					if (Os[l][id_v]) {
						int diff_v = i_ - Pos[l][id_v];
						if (diff_u < diff_v) {
							return 6;
						}
					}
				}
			}
		}

		if (!range_feasibility_original_node(I, S, l, i, bc))
			return 7;

		return 0;
	}

	if (range_feasibility_incremental_node(I, S, l, i, bc))
		return 0;

	return 1;

}

unsigned GRASP_v1::find_feasible_position_node(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned bc) {

	if (I.isOriginalNode(i, l))

		return find_feasible_position_original_node(I, S, l, i, bc);

	return find_feasible_position_incremental_node(I, S, l, i, bc);

}

unsigned GRASP_v1::find_feasible_position_original_node(HDAG &I, HDAG &S,
		unsigned l, unsigned i, unsigned bc) {

	auto k = getk();
	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();

	unsigned id_u = IDs[l][i];
	unsigned origin_pos = Pos[l][id_u];
	unsigned last_elem = LEVELS[l].size() - 1;

	if (!is_position_free(S, l, bc))
		bc = find_free_position(I, S, l, i, bc);

	unsigned violation = check_violation(I, S, l, i, bc);

	unsigned i_ = bc;

	while (true) {

		switch (violation) {
		case (0):
			return i_;
			break;
		case (1):
			// bc < p(u)
			i_ = origin_pos;
			while (!is_position_free(S, l, i_))
				i_++;

			violation = check_violation(I, S, l, i, i_);
			break;
		case (2):
			// bc > p(u) + k
			i_ = origin_pos + k;
			while (!is_position_free(S, l, i_))
				i_--;

			violation = check_violation(I, S, l, i, i_);

			break;
		case (3):
			// distance from p'(u) and last position does not consent
			// all remaining original nodes v such that p(v) > p(u)

			i_--;
			while (!is_position_free(S, l, i_))
				i_--;

			violation = check_violation(I, S, l, i, i_);
			break;
		case (4):
			// Exists v | p'(v) < bc and p(v) > p(u)
			i_--;
			while (!is_position_free(S, l, i_))
				i_--;

			violation = check_violation(I, S, l, i, i_);
			break;
		case (5):
			// Exists v | p'(v) > bc and p(v) < p(u)

			i_++;
			while (!is_position_free(S, l, i_))
				i_++;

			violation = check_violation(I, S, l, i, i_);
			break;
		case (6):
			// Exists v | p'(v) < bc and p'(v) - p(v) > p'(u) - p(u)
			i_++;
			while (!is_position_free(S, l, i_))
				i_++;
			violation = check_violation(I, S, l, i, i_);
			break;
		case (7):
			unsigned lower_limit = max(0, int((int)origin_pos - (int)k));
			unsigned upper_limit = min((int)last_elem, int((int)origin_pos + (int)k));
			unsigned initial_val = i_;
			if (abs((int)i_ - (int)lower_limit) < abs((int)i_ - (int)upper_limit)) {
				for (unsigned i_ = lower_limit; i_ <= upper_limit; ++i_) {
					if (i_ != initial_val) {
						if (i_ >= 0 && i_ < LEVELS[l].size()) {
							if (is_position_free(S, l, i_)) {
								violation = check_violation(I, S, l, i,
										i_);
								if (violation == 0)
									return i_;
							}
						}
					}
				}
			} else {
				for (unsigned i_ = upper_limit; i_ >= lower_limit; --i_) {
					if (i_ != initial_val) {
						if (i_ >= 0 && i_ < LEVELS[l].size()) {
							if (is_position_free(S, l, i_)) {
								violation = check_violation(I, S, l, i,
										i_);
								if (violation == 0)
									return i_;
							}
						}
					}
				}
			}
		}

	}
	return i_;

}

unsigned GRASP_v1::find_feasible_position_incremental_node(HDAG &I, HDAG &S,
		unsigned l, unsigned i, unsigned bc) {

	auto &LEVELS = I.getLEVELS();

	unsigned new_bc = bc;
	bool direction = true;
	unsigned up = 1;
	unsigned down = 1;

	if (!is_position_free(S, l, bc))
		new_bc = find_free_position(I, S, l, i, bc);

	unsigned violation = check_violation(I, S, l, i, new_bc);

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
			if (is_position_free(S, l, new_bc)) {

				violation = check_violation(I, S, l, i, new_bc);

			}
		}

	}

	return new_bc;

}

void GRASP_v1::add_first_node_to_solution(HDAG &I, HDAG &S, unsigned l,
		unsigned i, MTRand &rng) {

	auto k = getk();
	auto &LEVELS = I.getLEVELS();
	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &Os = I.getOs();

	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	auto &CL = getCL();

	unsigned pos = LEVELS[l].size();

	// If the node is an original node:
	if (I.isOriginalNode(i, l)) {

		pos = rng.randInt(k) + Pos[l][IDs[l][i]];

		unsigned difference_by_upmost = upmost_original[l] - Pos[l][IDs[l][i]];

		int upper_limit = (LEVELS[l].size() - 1) - difference_by_upmost;

		while ((int) pos > upper_limit) {

			pos = rng.randInt(k) + Pos[l][IDs[l][i]];
			upper_limit = (LEVELS[l].size() - 1) - pos;

		}

		original_to_assign[l]--;

	} // Otherwise
	else {

		// Special case
		if (k == 0) {

			unsigned i_;
			for (i_ = 0; i_ < LEVELS[l].size(); ++i_) {
				if (!I.isOriginalNode(i_, l))
					break;
			}
			pos = rng.randInt((LEVELS[l].size() - 1) - (i_ + 1)) + (i_ + 1);

		} // Insert the node in a random position in the level
		else {

			pos = rng.randInt(LEVELS[l].size() - 1);

		}

	}

	// The element is added in solution at the positon
	// specified by "pos"
	unsigned id_u = IDs[l][i];
	IDs_S[l][pos] = id_u;
	Pos_S[l][id_u] = pos;
	Os_S[l][pos] = Os[l][i];
	if (I.isOriginalNode(i, l)) {
		original_to_assign[l]--;
	}

	// Remove from CL
	CL[l][i] = 0;

	gettoAssign()[l]--;
	getTotalToAssign()--;

	}

bool GRASP_v1::range_feasibility_original_node(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned bc) {

	auto &LEVELS = I.getLEVELS();
	auto &Pos = I.getPos();
	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	// ---> Find first down position occupied by an original node
	int free_position_up = 0;
	int free_position_down = 0;
	int original_to_insert_up = 0;
	int original_to_insert_down = 0;
	int lower, upper;
	bool lower_node_found = false;
	bool upper_node_found = false;
	for (lower = bc - 1; lower >= 0; --lower) {
		if (IDs_S[l][lower] != LEVELS_S[l].size()) {
			if (Os_S[l][lower]) {
				lower_node_found = true;
				break;
			}
		}
		else
			free_position_down++;
	}
	// ---> Find first up position occupied by an original node
	for (upper = bc + 1; upper < (int)LEVELS[l].size(); ++upper) {
		if (IDs_S[l][upper] != LEVELS_S[l].size()) {
			if (Os_S[l][upper]) {
				upper_node_found = true;
				break;
			}
		}
		else
			free_position_up++;
	}

	int id_lower = LEVELS[l].size();
	int id_upper = LEVELS[l].size();
	int original_position_lower;
	int original_position_upper;
	if (lower_node_found) {
		id_lower = IDs_S[l][lower];
		original_position_lower = Pos[l][id_lower] + 1;
	}
	else {
		original_position_lower = 0;
	}
	if (upper_node_found) {
		id_upper = IDs_S[l][upper];
		original_position_upper = Pos[l][id_upper] - 1;
	}
	else {
		original_position_upper = LEVELS[l].size() - 1;
	}

	unsigned id = IDs[l][i];
	unsigned pos_id = Pos[l][id];

	// In the considered range there are also a lower original node
	// and an upper original node
	if (id_lower != (int)LEVELS[l].size() && id_upper != (int)LEVELS[l].size()) {

		for (int i_ = original_position_lower; i_ < (int)pos_id; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_down++;
			}
		}
		if (free_position_down < original_to_insert_down) return false;

		for (int i_ = (int)pos_id + 1; i_ <= original_position_upper; ++i_) {
			if (Os[l][i_])
				original_to_insert_up++;
		}
		if (free_position_up < original_to_insert_up) return false;
	}
	// In the considered range there is only a lower original node
	else if (id_lower != (int)LEVELS[l].size() && id_upper == (int)LEVELS[l].size()) {

		for (int i_ = original_position_lower; i_ < (int)pos_id; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_down++;
			}
		}
		if (free_position_down < original_to_insert_down) return false;

		for (int i_ = (int)pos_id + 1; i_ < original_position_upper; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_up++;
			}
		}
		if (free_position_up < original_to_insert_up) return false;

		unsigned first_free_pos = bc + 1;
		for (unsigned i__ = (int)pos_id + 1;  i__ <= (unsigned)original_position_upper; ++i__) {
			if (Os[l][i__]) {
				while (first_free_pos < LEVELS[l].size() && !is_position_free(S, l, first_free_pos))
				{ first_free_pos++; }
				if (first_free_pos > i__ + getk()) {
					return false;
				}
				first_free_pos++;
			}
		}

	}
	// In the considered range there is only an upper original node
	else if (id_lower == (int)LEVELS[l].size() && id_upper != (int)LEVELS[l].size()) {

		for (int i_ = 0; i_ < (int)pos_id; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_down++;
			}
		}
		if (free_position_down < original_to_insert_down) return false;

		for (int i_ = (int)pos_id + 1; i_ <= original_position_upper; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_up++;
			}
		}
		if (free_position_up < original_to_insert_up) return false;

	}
	// In the considered range both lower and upper original node don't exist
	else {

		for (int i_ = 0; i_ < (int)pos_id; ++i_) {
			if (Os[l][i_]) {
				original_to_insert_down++;
			}
		}
		if (free_position_down < original_to_insert_down) return false;


		for (int i_= (int)pos_id + 1; i_ <= original_position_upper; ++i_){
			if (Os[l][i_]) {
				original_to_insert_up++;
			}
		}
		if (free_position_up < original_to_insert_up) return false;

		unsigned first_free_pos = bc + 1;
		for (unsigned i__ = (int)pos_id + 1;  i__ <= (unsigned)original_position_upper; ++i__) {
			if (Os[l][i__]) {
				while (first_free_pos < LEVELS[l].size() && !is_position_free(S, l, first_free_pos))
				{ first_free_pos++; }
				if (first_free_pos > i__ + getk()) {
					return false;
				}
				first_free_pos++;
			}
		}


	}

	if (free_position_down >= original_to_insert_down)
		if (free_position_up >= original_to_insert_up)
			return true;

	return false;

}

bool GRASP_v1::range_feasibility_incremental_node(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned bc) {

	auto k = getk();
	auto &LEVELS = I.getLEVELS();
	auto &Pos = I.getPos();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	// ---> Find first down position occupied by an original node
	int free_position = 1;
	int original_to_insert = 0;
	int lower, upper;
	bool lower_node_found = false;
	bool upper_node_found = false;
	bool lower_incremental_found = false;
	bool upper_incremental_found = false;
	int remaining_shift = k;
	int first_feasibile_position = 0;
	int free_position_down = 0;
	for (lower = bc - 1; lower >= 0; --lower) {
		if (IDs_S[l][lower] != LEVELS_S[l].size()) {
			if (Os_S[l][lower]) {
				lower_node_found = true;
				break;
			}
			else {
				if (!lower_incremental_found) {
					lower_incremental_found = true;
				}
				if (lower <= (int)upmost_original[l]) {
					remaining_shift--;
				}
				if (remaining_shift <= 0) {
					first_feasibile_position = upmost_original[l] + getk() + 1;
				}
			}
		} else
			free_position++;
	}
	free_position_down = free_position - 1;
	// ---> Find first up position occupied by an original node
	for (upper = bc + 1; upper < (int) LEVELS[l].size(); ++upper) {
		if (IDs_S[l][upper] != LEVELS_S[l].size()) {
			if (Os_S[l][upper]) {
				upper_node_found = true;
				break;
			}
			else {
				if (!upper_incremental_found) {
					upper_incremental_found = true;
				}
				if (upper <= (int)upmost_original[l]) {
					remaining_shift--;
				}
				if (remaining_shift <= 0) {
					first_feasibile_position = upmost_original[l] + getk() + 1;
				}
			}
		} else
			free_position++;
	}

	int id_lower = LEVELS[l].size();
	int id_upper = LEVELS[l].size();
	int original_position_lower;
	int original_position_upper;
	if (lower_node_found) {
		id_lower = IDs_S[l][lower];
		original_position_lower = Pos[l][id_lower] + 1;
	} else {
		original_position_lower = 0;
	}
	if (upper_node_found) {
		id_upper = IDs_S[l][upper];
		original_position_upper = Pos[l][id_upper] - 1;
	} else {
		original_position_upper = LEVELS[l].size() - 1;
	}

	// In the considered range there are both a lower original node
	// and an upper original one
	if (id_lower != (int) LEVELS[l].size()
			&& id_upper != (int) LEVELS[l].size()) {

		for (int i_ = original_position_lower; i_ <= original_position_upper;
				++i_) {
			if (Os[l][i_]) {
				original_to_insert++;
			}
		}

	}
	// In the considered range there is only a lower original node
	else if (id_lower != (int) LEVELS[l].size()
			&& id_upper == (int) LEVELS[l].size()) {

		unsigned first_free_pos = bc + 1;
		for (int i_ = original_position_lower; i_ < original_position_upper;
				++i_) {
			if (Os[l][i_]) {
				original_to_insert++;
				if (free_position_down > 0) {
					free_position_down--;
				} else if (free_position_down == 0) {
					unsigned id_i_ = I.getIDs()[l][i_];
					while(first_free_pos < LEVELS[l].size() &&
							!is_position_free(S, l, first_free_pos))
					{first_free_pos++;}
					if (first_free_pos > Pos[l][id_i_] + getk()) {
						return false;
					}
					first_free_pos++;
				}
			}
		}
	}
	// In the considered range there is only an upper original node
	else if (id_lower == (int) LEVELS[l].size()
			&& id_upper != (int) LEVELS[l].size()) {

		for (int i_ = 0; i_ <= original_position_upper; ++i_) {
			if (Os[l][i_]) {
				original_to_insert++;
			}
		}
	}
	// In the considered range both lower and upper original node don't exist
	else {

		if ((int)bc < first_feasibile_position) return false;

		// Check if exists a lower incremental node and an upper incremental node
		if (lower_incremental_found && upper_incremental_found) {
			// Exists both a lower incremental node and an upper incremental one

			unsigned max_feasibile_pos_upmost_original;
			if (upmost_original[l] + k < LEVELS[l].size() - 1) {
				max_feasibile_pos_upmost_original = upmost_original[l] + k;
			}
			else {
				max_feasibile_pos_upmost_original = LEVELS[l].size() - 1;
			}
			free_position_down = 0;

			for (unsigned i_ = max_feasibile_pos_upmost_original; i_ > 0; --i_) {

				if (is_position_free(S, l, i_) && i_ != bc)
					free_position_down++;

			}

			if (bc != 0 && is_position_free(S, l, 0)) free_position_down++;

			if (free_position_down > (int)(upmost_original[l]))
				return true;
			return false;

		}
		else if (lower_incremental_found && !upper_incremental_found) {
			// Exists only a lower incremental node

			unsigned max_feasibile_pos_upmost_original;
			if (upmost_original[l] + k < LEVELS[l].size() - 1) {
				max_feasibile_pos_upmost_original = upmost_original[l] + k;
			}
			else {
				max_feasibile_pos_upmost_original = LEVELS[l].size() - 1;
			}
			free_position_down = 0;

			for (unsigned i_ = max_feasibile_pos_upmost_original; i_ > 0; --i_) {

				if (is_position_free(S, l, i_) && i_ != bc)
					free_position_down++;

			}

			if (bc != 0 && is_position_free(S, l, 0)) free_position_down++;

			if (free_position_down > (int)(upmost_original[l]))
				return true;
			return false;

		}
		else if (!lower_incremental_found && upper_incremental_found) {
			// Exists only an upper incremental node

			unsigned max_feasibile_pos_upmost_original;
			if (upmost_original[l] + k < LEVELS[l].size() - 1) {
				max_feasibile_pos_upmost_original = upmost_original[l] + k;
			}
			else {
				max_feasibile_pos_upmost_original = LEVELS[l].size() - 1;
			}
			free_position_down = 0;

			for (unsigned i_ = max_feasibile_pos_upmost_original; i_ > 0; i_--) {

				if (is_position_free(S, l, i_) && i_ != bc)
					free_position_down++;

			}

			if (bc != 0 && is_position_free(S, l, 0)) free_position_down++;


			if (free_position_down > (int)(upmost_original[l]))
				return true;
			return false;

		}
		else {
			// Both a lower incremental node and an upper incremental one don't exist
			if (k > 0) return true;
			else {
				if (bc <= upmost_original[l]) {
					return false;
				}
				return true;
			}
		}
	}

	if (free_position > original_to_insert)
		return true;

	return false;

}

unsigned GRASP_v1::find_free_position(HDAG &I, HDAG &S, unsigned l, unsigned i,
		unsigned bc) {

	unsigned new_bc;
	bool direction = true;
	unsigned up = 1;
	unsigned down = 1;
	auto &LEVELS = I.getLEVELS();

	while (true) {

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

		if (new_bc >= 0 && new_bc < LEVELS[l].size())
			if (is_position_free(S, l, new_bc))
				break;

	}

	return new_bc;

}
