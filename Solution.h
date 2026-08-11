/*
 * Solution.h
 *
 *  Created on: 22 mar 2017
 *      Author: antonio
 */

#include "HDAG.h"

#ifndef SOLUTION_H_
#define SOLUTION_H_

class Solution {
public:

	virtual ~Solution() {};

	virtual std::vector<unsigned> compute_level_gap_vector(HDAG &G) = 0;
	virtual std::pair<unsigned, unsigned> rescale_to_original_and_get_position(
				HDAG &G, std::vector<unsigned> lev_gap, unsigned u) = 0;

	virtual double getCost(HDAG &G) = 0;

	virtual void build_edge_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap) = 0;

	virtual void build_node_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap) = 0;

	virtual bool edge_pair_already_selected(unsigned i, unsigned j, unsigned k, unsigned l) = 0;

	virtual void add_edge_pair(unsigned i, unsigned j, unsigned k, unsigned l) = 0;

	virtual bool node_pair_already_selected(unsigned i, unsigned j) = 0;

	virtual void add_node_pair(unsigned i, unsigned j) = 0;

	virtual unsigned find_node_pair(unsigned i, unsigned j) = 0;

	virtual std::vector<std::pair<unsigned, unsigned> > &getN() { return N; }

	virtual std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > &getM() { return M; }


private:

	// Structure of pair of nodes
	std::vector<std::pair<unsigned, unsigned> > N;

	// Structure of pair of edges
	std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > M;

};

#endif /* SOLUTION_H_ */
