/*
 * SolutionIMLCM.h
 *
 *  Created on: 03 apr 2017
 *      Author: antonio
 */

#include "Solution.h"

#ifndef SOLUTIONIMLCM_H_
#define SOLUTIONIMLCM_H_

using namespace std;

class SolutionIMLCM : public Solution{
public:
	SolutionIMLCM();
	virtual ~SolutionIMLCM();

	std::vector<unsigned> compute_level_gap_vector(HDAG &G);
	std::pair<unsigned, unsigned> rescale_to_original_and_get_position(
				HDAG &G, std::vector<unsigned> lev_gap, unsigned u);
	double getCost(HDAG &G);

	double getCost_afterSwap(HDAG &G, unsigned l, double cost);

	double getCost_beforeOneSwap(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost);
	double getCost_afterOneSwap(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost);
	double getCost_beforeOneShift(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost);
	double getCost_afterOneShift(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned i_, double cost);

	void build_edge_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_node_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	bool edge_pair_already_selected(unsigned i, unsigned j, unsigned k, unsigned l);

	void add_edge_pair(unsigned i, unsigned j, unsigned k, unsigned l);

	bool node_pair_already_selected(unsigned i, unsigned j);

	void add_node_pair(unsigned i, unsigned j);

	unsigned find_node_pair(unsigned i, unsigned j);

	void initialize_costsForLevels(unsigned l) {
		costsForLevels.resize(l);
		bestsCostsForLevels.resize(l);
	};

	vector<double> &get_costsForLevels() { return costsForLevels; };
	vector<double> &get_bestsCostsForLevels() { return bestsCostsForLevels; };

private:

	vector<vector<pair<unsigned, unsigned> > > NODES;

	vector<vector<pair<pair<unsigned, unsigned>, pair<unsigned, unsigned> > > > ARCS;

	vector<double> costsForLevels;

	vector<double> bestsCostsForLevels;

};

#endif /* SOLUTIONIMLCM_H_ */
