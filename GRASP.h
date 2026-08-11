/*
 * GRASP.h
 *
 *  Created on: 23 mag 2017
 *      Author: antonio
 */

#include <vector>
#include "HDAG.h"
#include "lib/MTRand.h"
#include "SolutionIMLCM.h"
#include "chrono"
#include <string.h>
#include "PathRelinking.h"

using namespace std;

#ifndef GRASP_H_
#define GRASP_H_

class GRASP {
public:
	GRASP();
	virtual ~GRASP();

	virtual void algorithm(HDAG &I, int argc, char *argv[]) = 0;

	virtual void run(HDAG &I, HDAG &S) = 0;

	virtual void construction(HDAG &I, HDAG &S) = 0;

	virtual void allocate_and_build_all_structure(HDAG &I);
	virtual void initialize_all_iteration_structure(HDAG &I);
	virtual void delete_all_iterative_structure(HDAG &I);

	virtual std::vector<std::vector<unsigned>> &getCL() {
		return CL;
	};
	virtual std::vector<unsigned> &gettoAssign() {
		return toAssign;
	};
	std::vector<std::vector<std::pair<unsigned, unsigned>>>&getRCL() {return RCL;};

	void setk(unsigned x) {k = x;};
	void setalpha(double x) {alpha = x;};
	void setseed(unsigned x) {seed = x;};
	void setmaxIter(unsigned x) {max_iter = x;};
	void setTimeLimit(double t) {time_limit = t;};
	void setlocalsearch(char *x) { ls = x; };
	void setTotalToAssign(unsigned x) {total_to_assign = x;};
	void setTotalInRCL (unsigned x) {total_in_RCL = x;};

	const double getalpha() {return alpha;};
	const unsigned getseed() {return seed;};
	const unsigned getmaxIter() {return max_iter;}
	const char* getlocalsearch() {return ls;}
	const double getTimeLimit() {return time_limit;};
	unsigned &getTotalToAssign() {return total_to_assign;};
	unsigned &getTotalInRCL() {return total_in_RCL;};
	unsigned getk() {return k;};
	unsigned getrhomax() {return rhomax;};
	unsigned &getmaxdeg() {return maxdeg;};
	vector<vector<unsigned>> &getRhoCL() {return RhoCL;};
	vector<unsigned> &getresidual_shift() {return residual_shift;};
	vector<unsigned> &getupmost_inserted_node() {return upmost_inserted_node;};
	vector<vector<unsigned>> &getOD() {return OD;};
	vector<std::pair<unsigned, unsigned>> &getv_maxoutdeg() {return v_maxoutdeg;};

	void compute_degree_matrix(HDAG &I);

	void build_maxdegree_node_structure(HDAG &I, HDAG &S);

	pair<unsigned, unsigned> extract_maxdegree_node(HDAG &instance, MTRand &rng);

	void add_link(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc);

	void update_rho(HDAG &I, HDAG &S, unsigned l, unsigned i);

	void update_rhomax(HDAG &I, HDAG &S, unsigned l, unsigned i);

	void update_RCL(HDAG &I, HDAG &S, unsigned l, unsigned i);

	bool is_already_in_RCL(unsigned l, unsigned i);

	void build_RCL(HDAG &instance, HDAG &solution);

	void add_in_RCL(unsigned l, unsigned i) {
		pair<unsigned, unsigned> tmp;
		tmp.first = l;
		tmp.second = i;
		RCL[l].push_back(tmp);
		total_in_RCL++;
	}

	unsigned compute_potential_degree(HDAG &instance, HDAG &solution,
			unsigned l, unsigned i);

	bool is_in_CL(HDAG &instance, unsigned l, unsigned id);

	pair<unsigned, unsigned> extract_from_RCL(MTRand &rng);

	bool is_position_free(HDAG &solution, unsigned l, unsigned bc);

	void compute_rhomax(HDAG &I, HDAG &S);

	// ABSTRACT METHOD for CONSTRUCTION PHASE
	virtual void allocate_solution(HDAG &I, HDAG &S) = 0;
	virtual void add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng) = 0;
	virtual unsigned check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) = 0;
	virtual unsigned find_feasible_position_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) = 0;
	virtual unsigned find_feasible_position_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) = 0;
	virtual unsigned find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) = 0;

	void localsearch(HDAG &I, HDAG &S);

	/**************************************** FOR TABU ***************************************/
	void build_TABU_POSITION_MATRIX(HDAG &I);
	vector<vector<vector<unsigned>>> &get_TABU_POSITION_MATRIX()
			{ return TABU_POSITION_MATRIX; }
	void update_TABU_POSITION_MATRIX(HDAG &I, HDAG &S, unsigned iter);
	void localsearch_TABU(HDAG &I, HDAG &S, unsigned iter);
	bool is_TABU_POSITION(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned pos, unsigned iter) {

		if (iter - TABU_POSITION_MATRIX[l][i][pos] < getk()) {
			return true;
		}

		return false;
	}

	void compute_tenure(HDAG &I) {

		auto &INCREMENTAL = I.getINCREMENTAL();
		int min_incremental_number = INT_MAX;
		for (unsigned l = 0; l < INCREMENTAL.size(); ++l) {
			if (INCREMENTAL[l].size() < (unsigned)min_incremental_number)
				min_incremental_number = INCREMENTAL[l].size();
		}

		tenure_vertex_position = floor((double)(0.75 * (double)min_incremental_number));

	}

	void build_TABU_VERTEX_MATRIX(HDAG &I);
	vector<vector<unsigned>> &get_TABU_VERTEX_MATRIX()
			{ return TABU_VERTEX_MATRIX; }
	bool is_TABU_VERTEX (unsigned l, unsigned i, unsigned iter) {

		if (iter - TABU_VERTEX_MATRIX[l][i] < (unsigned)tenure_vertex_position) {

			return true;

		}

		return false;

	}
	/******************************************************************************************/

	void allocate_best_solution(HDAG &I);
	void copy_best_solution(HDAG &S);

	HDAG &get_best_solution() {return best_solution;}
	double &get_best_cost() {return best_cost;}
	double &get_time_to_best() {return time_to_best;}

	void set_best_cost(double x) {best_cost = x;}
	void set_time_to_best(double x) {time_to_best = x;}

	chrono::time_point<chrono::system_clock> &get_start_iteration() {return start_iteration;}
	chrono::time_point<chrono::system_clock> &get_end_iteration() {return end_iteration;}
	chrono::time_point<chrono::system_clock> &get_start_algorithm() {return start_algorithm;}
	double get_total_time_algorithm() {
		int elapsed =
				chrono::duration_cast<chrono::milliseconds>(end_iteration - start_algorithm).count();
		return (double)((double)elapsed / (double)1000);
	}

	void set_start_iteration(chrono::time_point<chrono::system_clock> x) {start_iteration = x;}
	void set_end_iteration(chrono::time_point<chrono::system_clock> x) {end_iteration = x;}
	void set_start_algorithm(chrono::time_point<chrono::system_clock> x) {start_algorithm = x;}

	void set_best_out_file(char *filename) { out_f_best = fopen(filename, "a"); };
	void set_complete_out_file(char *filename) { out_f_complete = fopen(filename, "a"); };
    void print_best_data_out_file() {
    	fprintf(out_f_best, "%d %0.2f %0.3f %0.3f %0.3f\n",
    			k, alpha, best_cost, get_total_time_algorithm(), time_to_best);
    }
    void print_complete_data_out_file(double actual_z, double current_t) {
    	fprintf(out_f_complete, "%0.2f %d %0.3f %0.3f\n",
    			alpha, seed, actual_z, current_t);
    }
	void close_out_file() { fclose(out_f_best); fclose(out_f_complete); };

	Elite_Set *getEliteSet() { return ES; };
	Path_Relinking *getPR() { return PR; };


private:

	// Separation value
	unsigned k;

	// Maximum degree indeg+outdeg between all the nodes in CL
	unsigned rhomax;

	// Maximum degree indeg+outdeg between all the nodes in the instance
	unsigned maxdeg;

	// Minimum and Maximum degree indeg+outdeg between all the nodes in the partial solution
	//unsigned minimum_deg_partial;
	//unsigned maximum_deg_partial;

	// is a value given by minimum_deg_partial + alpha(maximum_deg_partial - minimum_deg_partial)
	//double tau;

	// a matrix which specifies for each nodes in the partial solution the actual value
	// of its in+out degree
	//vector<vector<unsigned>> deg_partial;

	// alpha value for the constructive phase
	double alpha;

	// seed value for the generation of random value
	int seed;

	// maximum iterations value
	unsigned max_iter;

	// tenure for the TABU search based on the
	// vertex selection
	int tenure_vertex_position;

	// flag that indicates the type of local search selected
	char *ls;

	// value that establishes the number of nodes in the RCL
	unsigned total_in_RCL;

	// value that establishes the total number of node which have to be inserted
	// in the solution to obtain a complete solution for the problem
	unsigned total_to_assign;

	// file to store parameters about the best solution found
	FILE *out_f_best;

	// file to store all the parameters relative to all the iterations done
	FILE *out_f_complete;

	// vector that specifies the node that are not been assigned yet
	vector<unsigned> toAssign;

	// a matrix which specifies for each nodes in the instances its out-degree
	vector<vector<unsigned>> OD;

	// vector which contains all the nodes with maximum out degree to respect
	// the instance in input
	vector<std::pair<unsigned, unsigned>> v_maxoutdeg;

	// matrix which specifies the nodes which are in the candidate list yet
	vector<vector<unsigned>> CL;

	// matrix which specifies the nodes which are in the restricted candidate list
	vector<vector<pair<unsigned, unsigned>>> RCL;

	// matrix wich specifies the the degree in+out for each nodes in the RCL
	std::vector<std::vector<unsigned>> RhoCL;


	std::vector<unsigned> residual_shift;
	std::vector<unsigned> upmost_inserted_node;

	void phase_1(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost);
	void phase_1_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost, unsigned iter);

	void move(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost);
	void move_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost, unsigned iter);
	void move_TABU_2(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost, unsigned iter);

	unsigned find_best_shift_up(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve);
	unsigned find_best_shift_up_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);
	unsigned find_best_shift_up_TABU_2(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, double &cost_tabu, SolutionIMLCM &s);

	unsigned find_first_shift_up(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve);
	unsigned find_first_shift_up_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);

	void shift_down(HDAG &S, unsigned up_limit, unsigned dw_limit, unsigned l);

	unsigned find_best_shift_down(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve);
	unsigned find_best_shift_down_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);
	unsigned find_best_shift_down_TABU_2(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, double &cost_tabu, SolutionIMLCM &s);

	unsigned find_first_shift_down(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve);
	unsigned find_first_shift_down_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);

	void move_to_restablish_TABU(HDAG &I, HDAG &S, SolutionIMLCM &s, double &best_cost, unsigned iter);
	unsigned find_shift_up_to_restablish_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);
	unsigned find_shift_down_to_restablish_TABU(HDAG &I, HDAG &S, unsigned init, unsigned l,
			double &best_cost, SolutionIMLCM &s, bool &improve, unsigned iter);

	void shift_up(HDAG &S, unsigned dw_limit, unsigned up_limit, unsigned l);

	vector<vector<unsigned>> TABU_VERTEX_MATRIX;
	vector<vector<vector<unsigned>>> TABU_POSITION_MATRIX;

	HDAG best_solution;
	double best_cost;
	double time_to_best;
	double time_limit;

	Elite_Set *ES;
	Path_Relinking *PR;

	chrono::time_point<chrono::system_clock> start_algorithm;
	chrono::time_point<chrono::system_clock> start_iteration;
	chrono::time_point<chrono::system_clock> end_iteration;

};

#endif /* GRASP_H_ */
