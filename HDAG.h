/*
 * HDAG.h
 *
 *  Created on: 14 mar 2017
 *      Author: antonio
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <climits>

using namespace std;

#ifndef HDAG_H_
#define HDAG_H_

class HDAG {
public:
	HDAG();
	virtual ~HDAG();

	void read_instance(std::string filename);
	void read_jesusFormat(std::string filename);
	unsigned getLevNumber();
	unsigned getMaxDegree();

	std::vector< std::vector< std::vector< unsigned> > > &getLEVELS ();
	std::vector< std::vector< std::vector< unsigned> > > &getB_LEVELS ();
	std::vector< std::vector< unsigned > > &getINCREMENTAL();

	void setLEVELS(unsigned l);

	void allocateLEVELS (unsigned l);
	void allocateLevel (unsigned l, unsigned n);
	void allocateIDs (unsigned l);
	void allocateIDs (unsigned l, unsigned n);
	void allocateOs (unsigned l);
	void allocateOs (unsigned l, unsigned n);
	void allocatePos (unsigned l);
	void allocatePos (unsigned l, unsigned n);

	std::vector <std::vector< unsigned > > &getLevel(unsigned i);
	void setLevel(unsigned l, unsigned n);

	std::vector <std::vector< unsigned > > &getIDs();
	void setIDs(unsigned l, unsigned n);

	std::vector <std::vector< unsigned > > &getOs();
	void setOs(unsigned l, unsigned n);

	std::vector <std::vector< unsigned > > &getPos();
	void setPos(unsigned l, unsigned n);

	bool areCrossingEdge(unsigned i, unsigned j, unsigned i_, unsigned j_, unsigned l);
	bool isOriginalNode(unsigned i, unsigned l);
	bool isOriginalLink(unsigned i, unsigned j, unsigned l);
	bool isIncrementalNode(unsigned i, unsigned l);

	void swapPositions(unsigned p1, unsigned p2, unsigned l);

	unsigned compute_degree(unsigned l, unsigned i);

	void add_node_HDAG(HDAG &G, unsigned l, unsigned i, unsigned u);

	unsigned compute_barycenter(vector<vector<unsigned>> &CL, HDAG &I, unsigned l, unsigned i);

	void printHDAG(HDAG &G, std::string filename);

	void printHDAGgraphviz (HDAG &G, std::string filename);

	unsigned get_n_original_nodes(HDAG &G);

	int get_MIN_INCREMENTAL_NUMBER(HDAG &G) {
		int min = INT16_MAX;
		for (unsigned l = 0; l < this->N_INCREMENTAL.size(); ++l) {
			if ((int)N_INCREMENTAL[l] < min) {
				min = N_INCREMENTAL[l];
			}
		}
		return min;
	}

	unsigned get_total_nodes() { return total_nodes; }

	void set_total_nodes(unsigned x) { total_nodes = x; };

	void increase_total_nodes(unsigned x) { total_nodes+=x; };

	void copy_HDAG(HDAG &G);
	void copy_level(HDAG &G, unsigned l);
	void dealloc_HDAG();

	void set_cost(double x) { cost = x; };
	double get_cost() { return cost; };

	bool is_feasible(HDAG &I, unsigned k);

private:
	unsigned _l; // Number of levels in the hierarchical graph
	unsigned _n; // Max number of nodes for each levels;
	std::vector< std::vector< std::vector< unsigned> > > LEVELS; // levels of the HDAG
	std::vector< std::vector< std::vector< unsigned> > > B_LEVELS; // levels of the HDAG according the BS
	std::vector< std::vector< unsigned > > IDs; // identifiers of the vertices;
	std::vector< std::vector< unsigned > > Os; // index of original vertices 1=original - 0=not orignal;
	std::vector< std::vector< unsigned > > Pos; // position of each node;
	std::vector< std::vector< unsigned > > INCREMENTAL; // All the incremental node in the graphs;
	std::vector< unsigned > N_INCREMENTAL; // Number of incremental node for each levels
 	unsigned maxDegree;
 	unsigned total_nodes;
 	double cost;

};

#endif /* HDAG_H_ */
