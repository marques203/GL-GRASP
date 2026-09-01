/*
 * SolutionMLCM.h
 *
 *  Created on: 03 apr 2017
 *      Author: antonio
 */

#include "Solution.h"

#ifndef SOLUTIONMLCM_H_
#define SOLUTIONMLCM_H_

class SolutionMLCM : public Solution {
public:
	SolutionMLCM();
	virtual ~SolutionMLCM();


	std::vector<unsigned> compute_level_gap_vector(HDAG &G);
	std::pair<unsigned, unsigned> rescale_to_original_and_get_position(
				HDAG &G, std::vector<unsigned> lev_gap, unsigned u);
	double getCost(HDAG &G);


	void build_edge_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	void build_node_pair_structure(HDAG &G, std::vector<unsigned> &lev_gap);

	bool edge_pair_already_selected(unsigned i, unsigned j, unsigned k, unsigned l);

	void add_edge_pair(unsigned i, unsigned j, unsigned k, unsigned l);

	bool node_pair_already_selected(unsigned i, unsigned j);

	void add_node_pair(unsigned i, unsigned j);

	unsigned find_node_pair(unsigned i, unsigned j);

/*	std::vector<std::pair<unsigned, unsigned> > getN() { return N; }

	std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > getM() { return M; }

private:

	// Structure of pair of nodes
	std::vector<std::pair<unsigned, unsigned> > N;

	// Structure of pair of edges
	std::vector<std::pair<std::pair<unsigned, unsigned>, std::pair<unsigned, unsigned> > > M;
*/

};

#endif /* SOLUTIONMLCM_H_ */
