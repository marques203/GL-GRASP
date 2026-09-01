/*
 * SolutionMLCM.cpp
 *
 *  Created on: 03 apr 2017
 *      Author: antonio
 */

#include "SolutionMLCM.h"

using namespace std;

SolutionMLCM::SolutionMLCM() {
	// TODO Auto-generated constructor stub

}

SolutionMLCM::~SolutionMLCM() {
	// TODO Auto-generated destructor stub
}


std::vector<unsigned> SolutionMLCM::compute_level_gap_vector(HDAG &G) {

	const auto &LEVELS = G.getLEVELS();
	const auto &Os = G.getOs();

	vector<unsigned> v;

	v.push_back(0);

	unsigned gap = 0;

	for (unsigned l = 0; l < G.getLevNumber(); ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			if (Os[l][i])

				gap++;

		}

		v.push_back(gap);

	}

	return v;


}


double SolutionMLCM::getCost(HDAG &G) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();

	double cost = 0;

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

							if (G.isOriginalLink(u1, v1, l)
									&& G.isOriginalLink(u2, v2, l)) {

								if (G.areCrossingEdge(u1, v1, u2, v2, l)) {

									cost = cost + 1;

								}

							}

						}

					}

				}

			}

		}

	}

	return cost;

}


std::pair<unsigned, unsigned> SolutionMLCM::rescale_to_original_and_get_position(HDAG &G, std::vector<unsigned> lev_gap, unsigned u) {

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

void SolutionMLCM::build_edge_pair_structure(HDAG &G, vector<unsigned> &lev_gap) {

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

						if (G.isOriginalLink(u1, v1, l)
								&& G.isOriginalLink(u2, v2, l)) {

							u1 += lev_gap[l];
							v1 += lev_gap[l + 1];
							u2 += lev_gap[l];
							v2 += lev_gap[l + 1];

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

void SolutionMLCM::build_node_pair_structure(HDAG &G, vector<unsigned> &lev_gap) {

	const auto &LEVELS = G.getLEVELS();
	const auto &IDs = G.getIDs();

	for (unsigned l = 0; l < G.getLevNumber(); ++l) {

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {

			for (unsigned j = 0; j < LEVELS[l].size(); ++j) {

				for (unsigned k = 0; k < LEVELS[l].size(); ++k) {

					if (G.isOriginalNode(IDs[l][i], l)
							&& G.isOriginalNode(IDs[l][j], l)
							&& G.isOriginalNode(IDs[l][k], l)) {

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

}


bool SolutionMLCM::edge_pair_already_selected(unsigned i, unsigned j, unsigned k,
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

void SolutionMLCM::add_edge_pair(unsigned i, unsigned j, unsigned k, unsigned l) {

	std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > tmp;

	tmp.first.first = i;
	tmp.first.second = j;
	tmp.second.first = k;
	tmp.second.second = l;

	getM().push_back(tmp);

}

bool SolutionMLCM::node_pair_already_selected(unsigned i, unsigned j) {

	const auto &N = getN();

	for (unsigned t = 0; t < N.size(); ++t) {

		if (N[t].first == i && N[t].second == j)

			return true;

	}

	return false;

}

void SolutionMLCM::add_node_pair(unsigned i, unsigned j) {

	std::pair<unsigned, unsigned> tmp;
	tmp.first = i;
	tmp.second = j;

	getN().push_back(tmp);

}

unsigned SolutionMLCM::find_node_pair(unsigned i, unsigned j) {

	unsigned t;

	const auto &N = getN();

	for (t = 0; t < N.size(); ++t) {

		if (N[t].first == i && N[t].second == j)

			break;

	}

	return t;

}
