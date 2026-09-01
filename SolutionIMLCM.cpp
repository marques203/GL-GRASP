/*
 * SolutionIMLCM.cpp
 *
 *  Created on: 03 apr 2017
 *      Author: antonio
 */

#include "SolutionIMLCM.h"


using namespace std;

SolutionIMLCM::SolutionIMLCM() {
	// TODO Auto-generated constructor stub

}

SolutionIMLCM::~SolutionIMLCM() {
	// TODO Auto-generated destructor stub
}

std::vector<unsigned> SolutionIMLCM::compute_level_gap_vector(HDAG &G) {

	const auto &LEVELS = G.getLEVELS();

	vector<unsigned> v;

	v.push_back(0);

	unsigned gap = 0;

	for (unsigned l = 0; l < G.getLevNumber(); ++l) {

		gap += LEVELS[l].size();

		v.push_back(gap);

	}

	return v;

}

double SolutionIMLCM::getCost_afterSwap(HDAG &G, unsigned l, double cost) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();
	double new_cost = cost;

	if (l > 0) {

		l = l - 1;
 	    new_cost = new_cost - costsForLevels[l];
		costsForLevels[l] = 0;
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
							costsForLevels[l]++;
						}
					}
				}
			}
		}
		new_cost = new_cost + costsForLevels[l];
		l = l + 1;

	}

	if (l < LEVELS.size() - 1) {

	 	new_cost = new_cost - costsForLevels[l];
		costsForLevels[l] = 0;
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
							costsForLevels[l]++;
						}
					}
				}
			}
		}
		new_cost = new_cost + costsForLevels[l];
	}

	return new_cost;

}

double SolutionIMLCM::getCost_beforeOneSwap(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	unsigned x_min = min(i, i_);
	unsigned x_max = max(i, i_);

	if (l > 0) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			unsigned original_pos_x = Pos[l][id_x];
			for (unsigned j = 0; j < B_LEVELS[l][original_pos_x].size(); ++j) {
				unsigned id_j = B_LEVELS[l][original_pos_x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					unsigned original_pos_x_ = Pos[l][id_x_];
					for (unsigned j_ = 0;
							j_ < B_LEVELS[l][original_pos_x_].size(); ++j_) {
						unsigned id_j_ = B_LEVELS[l][original_pos_x_][j_];
						if (id_j != id_j_) {
							if (S.areCrossingEdge(id_j, id_x, id_j_, id_x_,
									l - 1)) {

								cost--;
							}
						}
					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			for (unsigned j = 0; j < LEVELS_S[l][x].size(); ++j) {
				// arc (x,j)
				unsigned id_j = LEVELS_S[l][x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					for (unsigned j_ = 0; j_ < LEVELS_S[l][x_].size(); ++j_) {
						// arc (x_, j_)
						unsigned id_j_ = LEVELS_S[l][x_][j_];
						if (id_j != id_j_) {

/*							cout << " (FS-minus) level = " << l <<
									" arcs = (" << id_x << ", " << id_j << ") - (" <<
									id_x_ << ", " << id_j_ << ")" << endl;*/

							if (S.areCrossingEdge(id_x, id_j, id_x_, id_j_, l)) {

/*								cout << " (FS-minus) level = " << l <<
										" arcs = (" << id_x << ", " << id_j << ") - (" <<
										id_x_ << ", " << id_j_ << ")" << endl;*/

								cost--;
							}
						}
					}
				}
			}
		}
	}

	return cost;

}

double SolutionIMLCM::getCost_afterOneSwap(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	unsigned x_min = min(i, i_);
	unsigned x_max = max(i, i_);

	if (l > 0) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			unsigned original_pos_x = Pos[l][id_x];
			for (unsigned j = 0; j < B_LEVELS[l][original_pos_x].size(); ++j) {
				unsigned id_j = B_LEVELS[l][original_pos_x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					unsigned original_pos_x_ = Pos[l][id_x_];
					for (unsigned j_ = 0;
							j_ < B_LEVELS[l][original_pos_x_].size(); ++j_) {
						unsigned id_j_ = B_LEVELS[l][original_pos_x_][j_];
						if (id_j != id_j_) {
							if (S.areCrossingEdge(id_j, id_x, id_j_, id_x_,
									l - 1)) {
								cost++;
							}
						}
					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {
		for (unsigned x = x_min; x <= x_max - 1; ++x) {
			unsigned id_x = IDs_S[l][x];
			for (unsigned j = 0; j < LEVELS_S[l][x].size(); ++j) {
				// arc (x,j)
				unsigned id_j = LEVELS_S[l][x][j];
				for (unsigned x_ = x + 1; x_ <= x_max; ++x_) {
					unsigned id_x_ = IDs_S[l][x_];
					for (unsigned j_ = 0; j_ < LEVELS_S[l][x_].size(); ++j_) {
						// arc (x_, j_)
						unsigned id_j_ = LEVELS_S[l][x_][j_];
						if (id_j != id_j_) {
							if (S.areCrossingEdge(id_x, id_j, id_x_, id_j_, l)) {

/*								cout << " (FS-plus) level = " << l <<
										" arcs = (" << id_x << ", " << id_j << ") - (" <<
										id_x_ << ", " << id_j_ << ")" << endl;*/

								cost++;
							}
						}
					}
				}
			}
		}
	}
	return cost;

}

double SolutionIMLCM::getCost_beforeOneShift(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	if (l > 0) {
		unsigned id_i = IDs_S[l][i];
		unsigned id_i_ = IDs_S[l][i_];
		unsigned original_pos_i = Pos[l][id_i];
		unsigned original_pos_i_ = Pos[l][id_i_];
		for (unsigned j = 0; j < B_LEVELS[l][original_pos_i].size(); ++j) {
			for (unsigned j_ = 0; j_ < B_LEVELS[l][original_pos_i_].size();
					++j_) {
				unsigned id_j = B_LEVELS[l][original_pos_i][j];
				unsigned id_j_ = B_LEVELS[l][original_pos_i_][j_];

				if (id_j != id_j_) {
					if (S.areCrossingEdge(id_j, id_i, id_j_, id_i_, l - 1)) {

						cost = cost - 1;

					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {

		for (unsigned j = 0; j < LEVELS_S[l][i].size(); ++j) {
			for (unsigned j_ = 0; j_ < LEVELS_S[l][i_].size(); ++j_) {
				// arc (i,j)
				unsigned id_i = IDs_S[l][i];
				unsigned id_j = LEVELS_S[l][i][j];

				// arc (i_,j_)
				unsigned id_i_ = IDs_S[l][i_];
				unsigned id_j_ = LEVELS_S[l][i_][j_];

				if (id_j != id_j_) {
					if (S.areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
						cost = cost - 1;

					}
				}
			}
		}
	}

	return cost;

}

double SolutionIMLCM::getCost_afterOneShift(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned i_, double cost) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();

	if (l > 0) {
		unsigned id_i = IDs_S[l][i];
		unsigned id_i_ = IDs_S[l][i_];
		unsigned original_pos_i = Pos[l][id_i];
		unsigned original_pos_i_ = Pos[l][id_i_];

		for (unsigned j = 0; j < B_LEVELS[l][original_pos_i].size(); ++j) {
			for (unsigned j_ = 0; j_ < B_LEVELS[l][original_pos_i_].size();
					++j_) {
				unsigned id_j = B_LEVELS[l][original_pos_i][j];
				unsigned id_j_ = B_LEVELS[l][original_pos_i_][j_];

				if (id_j != id_j_) {
					if (S.areCrossingEdge(id_j, id_i, id_j_, id_i_, l - 1)) {
						cost = cost + 1;

					}
				}
			}
		}
	}
	if (l < LEVELS.size() - 1) {

		for (unsigned j = 0; j < LEVELS_S[l][i].size(); ++j) {
			for (unsigned j_ = 0; j_ < LEVELS_S[l][i_].size(); ++j_) {
				// arc (i,j)
				unsigned id_i = IDs_S[l][i];
				unsigned id_j = LEVELS_S[l][i][j];

				// arc (i_,j_)
				unsigned id_i_ = IDs_S[l][i_];
				unsigned id_j_ = LEVELS_S[l][i_][j_];

				if (id_j != id_j_) {
					if (S.areCrossingEdge(id_i, id_j, id_i_, id_j_, l)) {
						cost = cost + 1;

					}
				}
			}
		}
	}

	return cost;

}

double SolutionIMLCM::getCost(HDAG &G) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();

	double cost = 0;

	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {

//		costsForLevels[l] = 0;
//		bestsCostsForLevels[l] = 0;

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

							cost++;

						}

					}

				}

			}

		}

	}

	/*
	for (unsigned l = 0; l < G.getLevNumber() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				unsigned level_gap = LEVELS[l].size();

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); i_++) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// Find edge (u1, v1)
						unsigned u1 = IDs[l][i];
						unsigned v1 = LEVELS[l][i][j];

						// Find edge (u2, v2)
						unsigned u2 = IDs[l][i_];
						unsigned v2 = LEVELS[l][i_][j_];

						if (u1 < (v1 + level_gap) && u2 < (v2 + level_gap)
								&& u1 < u2 && v1 != (v2 + level_gap)) {

							if (G.areCrossingEdge(u1, v1, u2, v2, l)) {

								cost = cost + 1;

							}

						}

					}

				}

			}

		}

	}
	*/

	return cost;

}


std::pair<unsigned, unsigned> SolutionIMLCM::rescale_to_original_and_get_position(HDAG &G, std::vector<unsigned> lev_gap, unsigned u) {

	pair<unsigned, unsigned> tmp;
	for (unsigned i = 0; i < lev_gap.size(); ++i) {
		if (u < lev_gap[i + 1]) {
			tmp.first = i;
			break;
		}
	}

	unsigned original_node = u - lev_gap[tmp.first];

	tmp.second = G.getPos()[tmp.first][original_node];

	return tmp;


}

void SolutionIMLCM::build_edge_pair_structure(HDAG &G, vector<unsigned> &lev_gap) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();

	// Select all pairs of edges
	for (unsigned l = 0; l < G.getLevNumber() - 1; ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				for (unsigned i_ = 0; i_ < LEVELS[l].size(); i_++) {

					for (unsigned j_ = 0; j_ < LEVELS[l][i_].size(); ++j_) {

						// Find edge (u1, v1)
						unsigned u1 = IDs[l][i];
						unsigned v1 = LEVELS[l][i][j];

						// Find edge (u2, v2)
						unsigned u2 = IDs[l][i_];
						unsigned v2 = LEVELS[l][i_][j_];

						u1 += lev_gap[l];
						v1 += lev_gap[l + 1];
						u2 += lev_gap[l];
						v2 += lev_gap[l + 1];

						if (u1 != u2 && v1 != v2) {

							if (!edge_pair_already_selected(u1, v1, u2, v2)) {

								if (v1 < v2) {

									if (!node_pair_already_selected(v1, v2))
										add_node_pair(v1, v2);


									if (!node_pair_already_selected(u1, u2))
										add_node_pair(u1, u2);

									add_edge_pair(u1, v1, u2, v2);

								} else if (v1 > v2) {

									if (!node_pair_already_selected(v2, v1))
										add_node_pair(v2, v1);

									if (!node_pair_already_selected(u1, u2))
										add_node_pair(u1, u2);

									add_edge_pair(u1, v1, u2, v2);

								}

							}

						}
					}

				}

			}

		}

	}

}

void SolutionIMLCM::build_node_pair_structure(HDAG &G, vector<unsigned> &lev_gap) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();

	for (unsigned l = 0; l < G.getLevNumber(); ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l].size(); ++j) {

				for (unsigned k = 0; k < LEVELS[l].size(); ++k) {

					unsigned u1 = IDs[l][i] + lev_gap[l];
					unsigned u2 = IDs[l][j] + lev_gap[l];
					unsigned u3 = IDs[l][k] + lev_gap[l];

					if (u1 < u2 && u2 < u3) {

						if (!node_pair_already_selected(u1, u2))
							add_node_pair(u1, u2);

						if (!node_pair_already_selected(u2, u3))
							add_node_pair(u2, u3);

						if (!node_pair_already_selected(u1, u3))
							add_node_pair(u1, u3);

					}

				}

			}

		}

	}

}


bool SolutionIMLCM::edge_pair_already_selected(unsigned i, unsigned j, unsigned k,
		unsigned l) {

	const auto &M = getM();

	for (unsigned t = 0; t < M.size(); ++t) {

		if ((M[t].first.first == i && M[t].first.second == j)
				&& (M[t].second.first == k && M[t].second.second == l))

			return true;

		if ((M[t].first.first == k && M[t].first.second == l)
				&& (M[t].second.first == i && M[t].second.second == j))

			return true;

	}

	return false;

}

void SolutionIMLCM::add_edge_pair(unsigned i, unsigned j, unsigned k, unsigned l) {

	std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > tmp;

	tmp.first.first = i;
	tmp.first.second = j;
	tmp.second.first = k;
	tmp.second.second = l;

	getM().push_back(tmp);

}

bool SolutionIMLCM::node_pair_already_selected(unsigned i, unsigned j) {

	const auto &N = getN();

	for (unsigned t = 0; t < N.size(); ++t) {

		if (N[t].first == i && N[t].second == j)

			return true;

	}

	return false;

}

void SolutionIMLCM::add_node_pair(unsigned i, unsigned j) {

	std::pair<unsigned, unsigned> tmp;
	tmp.first = i;
	tmp.second = j;

	getN().push_back(tmp);

}

unsigned SolutionIMLCM::find_node_pair(unsigned i, unsigned j) {

	unsigned t;

	const auto &N = getN();

	for (t = 0; t < N.size(); ++t) {

		if (N[t].first == i && N[t].second == j)

			break;

	}

	return t;

}
