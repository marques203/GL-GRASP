/*
 * HDAG.cpp
 *
 *  Created on: 14 mar 2017
 *      Author: antonio
 */

#include "HDAG.h"
#include "SolutionMLCM.h"
#include "SolutionIMLCM.h"
#include <math.h>
#include "string.h"

using namespace std;

HDAG::HDAG():_l(0), _n(0), maxDegree(0), total_nodes(0), cost(INT_FAST64_MAX) {
	// TODO Auto-generated constructor stub

}

void HDAG::read_instance(string filename) {

	ifstream f;
	f.open(filename.c_str());

	int n_levels = 0;
	int n_nodes = 0;
	int n_new_nodes = 0;
	int n_arcs = 0;

	try {

		string line;

		getline(f, line);

		line = line + "\n";

		std::istringstream iss(line);
		iss >> _l;

		n_levels = _l;

		LEVELS.resize(_l);
		B_LEVELS.resize(_l);
		IDs.resize(_l);
		Os.resize(_l);
		Pos.resize(_l);
		INCREMENTAL.resize(_l);
		N_INCREMENTAL.resize(_l);

		getline(f, line);

		line = line + "\n";

		iss.clear();
		iss.str(line);
		for (unsigned i = 0; i < _l; ++i) {
			unsigned x;
			iss >> x;
			LEVELS[i].resize(x);
			B_LEVELS[i].resize(x);
			IDs[i].resize(x);
			Os[i].resize(x);
			Pos[i].resize(x);
			total_nodes = total_nodes + x;
		}

		unsigned n_lev = 0;
		for (; n_lev < _l - 1; ++n_lev) {

			for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

				getline(f, line);

				line = line + "\n";

				iss.clear();
				iss.str(line);
				unsigned o, u, v;
				iss >> o >> u;
				IDs[n_lev][i] = u;
				Os[n_lev][i] = o;
				Pos[n_lev][u] = i;

				n_nodes++;

				if (o == 0) {
					INCREMENTAL[n_lev].push_back(u);
					N_INCREMENTAL[n_lev]++;
					n_new_nodes++;
				}

				while (1) {

					iss >> v;

					if (iss.eof() == 1)
						break;

					LEVELS[n_lev][i].push_back(v);

					n_arcs++;

				}

				unsigned deg = compute_degree(n_lev, i);
				if (deg > maxDegree)
					maxDegree = deg;

			}

		}

		for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

			getline(f, line);
			iss.clear();
			iss.str(line);
			unsigned o, u;
			iss >> o >> u;
			IDs[n_lev][i] = u;
			Os[n_lev][i] = o;
			Pos[n_lev][u] = i;

			if (o == 0) {
				INCREMENTAL[n_lev].push_back(u);
				N_INCREMENTAL[n_lev]++;
			}

		}

	} catch (ifstream::failure &e) {

		cerr << "Failure when reading configuration-file: " << e.what() << endl;
		exit (0);

	}

	// BUILD BACKWARD_STAR
	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned id = IDs[l][i];
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				unsigned id_ = LEVELS[l][i][j];
				unsigned pos_id_ = Pos[l + 1][id_];
				B_LEVELS[l + 1][pos_id_].push_back(id);
			}
		}
	}

}

void HDAG::read_jesusFormat(std::string filename) {

	ifstream f;
	f.open(filename.c_str());

	int n_levels = 0;
	int n_nodes = 0;
	int n_new_nodes = 0;
	int n_arcs = 0;

	try {

		string line;

		getline(f, line);

		line = line + "\n";

		std::istringstream iss(line);
		iss >> _l;

		n_levels = _l;

		LEVELS.resize(_l);
		B_LEVELS.resize(_l);
		IDs.resize(_l);
		Os.resize(_l);
		Pos.resize(_l);
		INCREMENTAL.resize(_l);
		N_INCREMENTAL.resize(_l);

		getline(f, line);

		line = line + "\n";

		iss.clear();
		iss.str(line);
		for (unsigned i = 0; i < _l; ++i) {
			unsigned x;
			iss >> x;
			LEVELS[i].resize(x);
			B_LEVELS[i].resize(x);
			IDs[i].resize(x);
			Os[i].resize(x);
			Pos[i].resize(x);
			total_nodes = total_nodes + x;
		}

		int position;
		int id_node;

		unsigned n_lev = 0;
		for (; n_lev < _l - 1; ++n_lev) {

			id_node = 0;

			for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

				getline(f, line);

				line = line + "\n";

				iss.clear();
				iss.str(line);
				unsigned o, u, v;
				iss >> o >> u;

				position = u;

				if (o == 0) {

					INCREMENTAL[n_lev].push_back(id_node);
					N_INCREMENTAL[n_lev]++;
					n_new_nodes++;

				}

				IDs[n_lev][position] = id_node;
				Os[n_lev][position] = o;
				Pos[n_lev][id_node] = position;

//				cout << " node " << id_node << " in position " << position << endl;
//				getchar();

				n_nodes++;
				id_node++;

				while (1) {

					iss >> v;

					if (iss.eof() == 1)
						break;

					LEVELS[n_lev][position].push_back(v);

					n_arcs++;

				}

				unsigned deg = compute_degree(n_lev, i);
				if (deg > maxDegree)
					maxDegree = deg;

			}

		}

		id_node = 0;

		for (unsigned i = 0; i < LEVELS[n_lev].size(); ++i) {

			getline(f, line);
			iss.clear();
			iss.str(line);
			unsigned o, u;
			iss >> o >> u;

			position = u;

			if (o == 0) {
				INCREMENTAL[n_lev].push_back(u);
				N_INCREMENTAL[n_lev]++;
			}

			IDs[n_lev][position] = id_node;
			Os[n_lev][position] = o;
			Pos[n_lev][id_node] = position;

			id_node++;

		}


	} catch (ifstream::failure &e) {

		cerr << "Failure when reading configuration-file: " << e.what() << endl;
		exit (0);

	}

	// BUILD BACKWARD_STAR
	for (unsigned l = 0; l < LEVELS.size() - 1; ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned id = IDs[l][i];
			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
				unsigned id_ = LEVELS[l][i][j];
				unsigned pos_id_ = Pos[l + 1][id_];
				B_LEVELS[l + 1][pos_id_].push_back(id);
			}
		}
	}

}

HDAG::~HDAG() {
	// TODO Auto-generated destructor stub
}

void HDAG::allocateLEVELS(unsigned l) {

	this->LEVELS.resize(l);
	this->_l = l;

}

void HDAG::allocateLevel(unsigned l, unsigned n) {

	this->LEVELS[l].resize(n);

}

void HDAG::allocateIDs(unsigned l) {

	this->IDs.resize(l);

}

void HDAG::allocateIDs(unsigned l, unsigned n) {

	this->IDs[l].resize(n, this->LEVELS[l].size());

}

void HDAG::allocateOs(unsigned l) {

	this->Os.resize(l);

}

void HDAG::allocateOs(unsigned l, unsigned n) {

	this->Os[l].resize(n, this->LEVELS[l].size());

}

void HDAG::allocatePos(unsigned l) {

	this->Pos.resize(l);

}

void HDAG::allocatePos(unsigned l, unsigned n) {

	this->Pos[l].resize(n, this->LEVELS[l].size());

}

std::vector< std::vector< std::vector< unsigned> > > &HDAG::getLEVELS() {

	return this->LEVELS;

}

std::vector< std::vector< std::vector< unsigned> > > &HDAG::getB_LEVELS() {

	return this->B_LEVELS;

}


unsigned HDAG::getLevNumber() {

	return this->LEVELS.size();

}

unsigned HDAG::getMaxDegree() {

	return this->maxDegree;

}

std::vector<std::vector<unsigned> > &HDAG::getLevel(unsigned i) {

	return this->LEVELS[i];

}

std::vector< std::vector< unsigned > > &HDAG::getINCREMENTAL() {

	return this->INCREMENTAL;

}

std::vector<std::vector<unsigned> > &HDAG::getIDs() {

	return this->IDs;

}

std::vector<std::vector<unsigned> > &HDAG::getOs() {

	return this->Os;

}

std::vector<std::vector<unsigned> > &HDAG::getPos() {

	return this->Pos;

}

bool HDAG::areCrossingEdge(unsigned i, unsigned j, unsigned i_, unsigned j_,
		unsigned l) {


	if ((this->Pos[l][i] < this->Pos[l][i_]) && (this->Pos[l + 1][j] > this->Pos[l + 1][j_]))
		return true;

	if ((this->Pos[l][i] > this->Pos[l][i_]) && (this->Pos[l + 1][j] < this->Pos[l + 1][j_]))
		return true;

	return false;

}

bool HDAG::isOriginalLink(unsigned i, unsigned j, unsigned l) {

	if (this->Os[l][this->Pos[l][i]] == 1 && this->Os[l + 1][this->Pos[l + 1][j]] == 1)
		return true;

	return false;

}

bool HDAG::isOriginalNode(unsigned i, unsigned l) {

	if (this->Os[l][i] == 1)
		return true;

	return false;

}

bool HDAG::isIncrementalNode(unsigned i, unsigned l) {

	if (!isOriginalNode(i, l))
		return true;

	return false;

}

void HDAG::swapPositions(unsigned p1, unsigned p2, unsigned l) {

	this->Pos[l][this->IDs[l][p1]] = p2;
	this->Pos[l][this->IDs[l][p2]] = p1;

	this->LEVELS[l][p1].swap(this->LEVELS[l][p2]);

	unsigned tmp = this->Os[l][p1];
	this->Os[l][p1] = this->Os[l][p2];
	this->Os[l][p2] = tmp;

	tmp = this->IDs[l][p1];
	this->IDs[l][p1] = this->IDs[l][p2];
	this->IDs[l][p2] = tmp;

	return;

}

unsigned HDAG::compute_degree(unsigned l, unsigned i) {

	auto &LEVELS = this->getLEVELS();
	unsigned deg = 0;

	if (l > 0) {

		for (unsigned i_ = 0; i_ < LEVELS[l - 1].size(); ++i_) {
			for (unsigned j_ = 0; j_ < LEVELS[l - 1][i_].size(); ++j_) {
				if 	(LEVELS[l - 1][i_][j_] == IDs[l][i]) {
					deg++;
					break;
				}
			}
		}
	}

	deg += LEVELS[l][i].size();

	return deg;

}

unsigned HDAG::compute_barycenter(vector<vector<unsigned>> &CL, HDAG &I, unsigned l, unsigned i) {

	auto &LEVELS = I.getLEVELS();
	auto &B_LEVELS = I.getB_LEVELS();
	auto &Pos = I.getPos();

	auto &Pos_S = this->getPos();

	double BS_bc = 0;
	double FS_bc = 0;
	unsigned BS_elem = 0;
	unsigned FS_elem = 0;

	// compute the barycenter in backward sense
	for (unsigned j = 0; j < B_LEVELS[l][i].size(); ++j) {
		unsigned id_j = B_LEVELS[l][i][j];
		unsigned original_pos_id_j = Pos[l - 1][id_j];
		// if the node is already in the new solution
		if (CL[l - 1][original_pos_id_j] == 0) {
			unsigned new_pos_id_j = Pos_S[l - 1][id_j];
			BS_bc = BS_bc + new_pos_id_j;
			BS_elem++;
		}
	}

	// compute the barycenter in forward sense
	for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {
		unsigned id_j = LEVELS[l][i][j];
		unsigned original_pos_id_j = Pos[l + 1][id_j];
		// if the node is already in the new solution
		if (CL[l + 1][original_pos_id_j] == 0) {
			unsigned new_pos_id_j = Pos_S[l + 1][id_j];
			FS_bc = FS_bc + new_pos_id_j;
			FS_elem++;
		}
	}

	if (BS_elem != 0)
		BS_bc = BS_bc / BS_elem;
	if (FS_elem != 0)
		FS_bc = FS_bc / FS_elem;

	if (BS_elem == 0 && FS_elem == 0) return 0;

	unsigned bc = (floor)((BS_bc + FS_bc) / 2);

	return bc;

/*	const auto &LEVELS = this->getLEVELS();
	const auto &IDs = this->getIDs();

	unsigned position_sum = 0;
	unsigned elems_summed = 0;

	double bc_forward = floor(LEVELS[l].size() / 2);
	double bc_backward = floor(LEVELS[l].size() / 2);

	// Compute the barycenter in forward sense
	if (l < LEVELS.size() - 1) {

		for (unsigned i_ = 0; i_ < LEVELS[l+1].size(); ++i_) {

			if (IDs[l+1][i_] != LEVELS[l+1].size()) {

				position_sum += i_;
				elems_summed++;

			}

		}

		if (elems_summed > 0)
			bc_forward = floor(position_sum / elems_summed);

	}

	position_sum = 0;
	elems_summed = 0;

	//Compute the baryceter in backward sense
	if (l > 0) {

		for (unsigned i_ = 0; i_ < LEVELS[l-1].size(); ++i_) {

			if (IDs[l-1][i_] != LEVELS[l-1].size()) {

				position_sum += i_;
				elems_summed++;

			}

		}

		if (elems_summed > 0)
			bc_backward = floor(position_sum / elems_summed);

	}

	unsigned bc = floor((bc_forward + bc_backward) / 2);

	return bc;*/

}

unsigned HDAG::get_n_original_nodes(HDAG &G) {

	auto &LEVELS = G.getLEVELS();
	unsigned n = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l)
		for (unsigned i = 0; i < LEVELS[l].size(); ++i)
			if (G.isOriginalNode(i, l))
				n++;

	return n;

}

void HDAG::printHDAG(HDAG &G, string filename) {

	ofstream f(filename.c_str());

	f << "\\documentclass{standalone}" << '\n';
	f << "\\usepackage{tikz}" << '\n';
	f << "\\usepackage{amsmath}" << '\n';
	f << "\\usepackage{amsfonts}" << '\n';
	f << "\\usepackage{comment}" << '\n';
	f << "\\begin{document}" << '\n';
	f << "\\begin{tikzpicture}[thick]" << '\n';
	f << "\\definecolor{gainsboro}{rgb}{0.86, 0.86, 0.86} " << "\n";

	// Print the nodes
	unsigned index_l = 0;

	unsigned max_dim_level = 0;
	double x_position_label_cost = ((G.getLevNumber() - 1) * 10) / 2;

	f << "% VERTICES" << '\n';

	for (unsigned i = 0; i < G.getLevNumber(); ++i) {

		int index_i = 0;

		f << "% Node of level " << i << '\n';

		if (G.getLevel(i).size() > max_dim_level)
			max_dim_level = G.getLevel(i).size();

		int node_in_layers = G.getLevel(i).size();

		for (unsigned j = 0; j < node_in_layers; ++j) {

			// revert print
			int x = (node_in_layers - 1) - j;

			if (G.Os[i][x] == 1) {

				f << "\\node [draw, circle, fill=gainsboro, minimum size = 0.7cm] (" << G.IDs[i][x] << "-" << i
						<< ") at (" << index_l << ", " << index_i << ") {\\sffamily\\small"
						<< G.IDs[i][x] + 1 << "};" << '\n';

			} else {

				if (G.IDs[i][j] == LEVELS[i].size()) {

					f << "\\node [draw, circle, fill=black, minimum size = 0.7cm, text=white] (" << G.IDs[i][x] << "-"
							<< i << ") at (" << index_l << ", " << index_i << ") {\\sffamily\\small x };" << '\n';

				}
				else {
					f << "\\node [draw, circle, fill=black, minimum size = 0.7cm, text=white] (" << G.IDs[i][x] << "-"
						<< i << ") at (" << index_l << ", " << index_i << ") {\\sffamily\\small"
						<< G.IDs[i][x] + 1 << "};" << '\n';
				}
			}

			index_i++;

		}

		f << '\n';

		index_l = index_l + 12;

	}

	f << "% EDGES" << '\n';

	//Print edges
	for (unsigned i = 0; i < G.getLevNumber() - 1; ++i) {

		f << "% Connection between levels " << i << " and " << i + 1 << '\n';

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			for (unsigned k = 0; k < G.getLevel(i)[j].size(); ++k) {

				// Find link (u, v)
				unsigned u = G.IDs[i][j];
				unsigned v = G.LEVELS[i][j][k];

				if (isOriginalLink(u, v, i)) {

					f << "\\draw[gainsboro] (" << G.IDs[i][j] << "-" << i << ") -- ("
							<< G.getLevel(i)[j][k] << "-" << i + 1 << ");"
							<< '\n';

				} else {

					f << "\\draw (" << G.IDs[i][j] << "-" << i
							<< ") -- (" << G.getLevel(i)[j][k] << "-" << i + 1
							<< ");" << '\n';

				}

			}

		}

	}
//	f << "\\node [font=\\Large] at (" << x_position_label_cost << ", " << max_dim_level + 1
//			<< ") {Total nodes = " << G.get_cost() << "};"
//			<< "'\n";

	f << "\\end{tikzpicture}" << '\n';
	f << "\\end{document}" << '\n';

	f.close();

	return;

}

void HDAG::printHDAGgraphviz (HDAG &G, std::string filename) {

	ofstream f(filename.c_str());

	f << "digraph{" << "\n";

	// Print the nodes
	unsigned index_l = 0;

	unsigned max_dim_level = 0;
	double x_position_label_cost = ((G.getLevNumber() - 1) * 10) / 2;


	for (unsigned i = 0; i < G.getLevNumber(); ++i) {

		int index_i = 0;


		if (G.getLevel(i).size() > max_dim_level)
			max_dim_level = G.getLevel(i).size();

		int node_in_layers = G.getLevel(i).size();

		for (unsigned j = 0; j < node_in_layers; ++j) {

			// revert print
			int x = (node_in_layers - 1) - j;

			if (G.Os[i][x] == 1) {

				f << "n" << G.IDs[i][x] << "a" << i << " [label=" << G.IDs[i][x] << ","
						<< "shape=circle,color=gray,style=filled,fontcolor=black,fontname=\"times-bold\" ];" << '\n';

			} else {

				if (G.IDs[i][j] == LEVELS[i].size()) {

					f << "n" << G.IDs[i][x] << "a" << i << " [label=" << G.IDs[i][x] << ","
							<< "shape=circle,color=black,style=filled,fontcolor=white,fontname=\"times-bold\"" << "];" << '\n';

				}
				else {

					f << "n" << G.IDs[i][x] << "a" << i << " [label=" <<  G.IDs[i][x] << ","
							<< "shape=circle,color=black,style=filled,fontcolor=white,fontname=\"times-bold\"" << "];" << '\n';

				}
			}

			index_i++;

		}

		f << '\n';

		index_l = index_l + 12;

	}

	//Print edges
	for (unsigned i = 0; i < G.getLevNumber() - 1; ++i) {

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			for (unsigned k = 0; k < G.getLevel(i)[j].size(); ++k) {

				// Find link (u, v)
				unsigned u = G.IDs[i][j];
				unsigned v = G.LEVELS[i][j][k];

				if (isOriginalLink(u, v, i)) {

//					f << "\\draw[gainsboro] (" << G.IDs[i][j] << "-" << i << ") -- ("
//							<< G.getLevel(i)[j][k] << "-" << i + 1 << ");"
//							<< '\n';

					f << "n" << G.IDs[i][j] << "a" << i << " -> " << "n" << G.getLevel(i)[j][k] << "a" << i + 1
							<< " [dir=none, color=\"gray\"]" << '\n';


				} else {

					f << "n" << G.IDs[i][j] << "a" << i << " -> " << "n" << G.getLevel(i)[j][k] << "a" << i + 1
							<< " [dir=none]" << '\n';

				}

			}

		}

	}

	f << "}" << '\n';

	f.close();

	return;


}

void HDAG::copy_HDAG(HDAG &G) {

	auto &LEVELS = G.getLEVELS();
	auto &IDs = G.getIDs();
	auto &Os = G.getOs();

	auto &IDs_new = this->getIDs();
	auto &Pos_new = this->getPos();
	auto &Os_new = this->getOs();

	unsigned ln, vn;

	ln = G.getLevNumber();

	this->allocateLEVELS(ln);
	this->allocateIDs(ln);
	this->allocateOs(ln);
	this->allocatePos(ln);

	for (unsigned l = 0; l < ln; ++l) {

		vn = LEVELS[l].size();

		this->allocateLevel(l, vn);
		this->allocateIDs(l, vn);
		this->allocateOs(l, vn);
		this->allocatePos(l, vn);

		for (unsigned i = 0; i < vn; ++i) {

			unsigned id_u = IDs[l][i];
			IDs_new[l][i] = id_u;
			Pos_new[l][id_u] = i;
			Os_new[l][i] = Os[l][i];

			for (unsigned j = 0; j < LEVELS[l][i].size(); ++j) {

				this->getLEVELS()[l][i].push_back(LEVELS[l][i][j]);

			}

		}

	}

	return;

}

void HDAG::copy_level(HDAG &G, unsigned l) {

	auto &LEVELS_target = G.getLEVELS();
	auto &IDs_target = G.getIDs();

	auto &IDs_source = this->getIDs();
	auto &Pos_source = this->getPos();

	for (unsigned i = 0; i < LEVELS_target[l].size(); ++i) {

		unsigned id_v_in_target = IDs_target[l][i];
		unsigned pos_v_in_target = i;

		unsigned id_v_in_source = id_v_in_target;
		unsigned pos_v_in_source = Pos_source[l][id_v_in_source];

		if (pos_v_in_target != pos_v_in_source) {

			this->swapPositions(pos_v_in_source, pos_v_in_target, l);

		}

	}

	return;

}

void HDAG::dealloc_HDAG() {

	for (unsigned l = 0; l < this->getLEVELS().size(); ++l) {
		for (unsigned i = 0; i < this->getLEVELS()[l].size(); ++i) {
			while (this->getLEVELS()[l][i].size() > 0)
				this->getLEVELS()[l][i].pop_back();
			delete &this->LEVELS[l][i];
		}
		delete &this->LEVELS[l];
	}
	delete this;

}

bool HDAG::is_feasible(HDAG &I, unsigned k) {

	const auto &LEVELS = I.getLEVELS();
	const auto &LEVELS_S = this->getLEVELS();

	const auto &IDs = I.getIDs();
	const auto &IDs_S = this->getIDs();

	const auto &Pos = I.getPos();
	const auto &Pos_S = this->getPos();

	const auto &Os = I.getOs();
	const auto &Os_S = this->getOs();

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			unsigned id = IDs[l][i];
			unsigned old_pos = i;
			unsigned new_pos = Pos_S[l][id];

			if ( (IDs_S[l][i] == LEVELS_S[l].size()) ||
				 (Pos_S[l][i] == LEVELS_S[l].size()) ||
				 ( Os_S[l][i] == LEVELS_S[l].size())    )
				return false;

			if (Os[l][i]) {
				if (new_pos > (old_pos + k))
					return false;
			}
			unsigned dim_FS_old = LEVELS[l][old_pos].size();
			unsigned dim_FS_new = LEVELS_S[l][new_pos].size();
			if (dim_FS_old != dim_FS_new)
				return false;
			bool neighborhood_found = false;
			for (unsigned j = 0; j < LEVELS[l][old_pos].size(); ++j) {
				unsigned neigh_old = LEVELS[l][old_pos][j];
				for (unsigned j_ = 0; j_ < LEVELS_S[l][new_pos].size(); ++j_) {
					if (neigh_old == LEVELS_S[l][new_pos][j_]) {
						neighborhood_found = true;
						break;
					}
				}
				if (!neighborhood_found)
					return false;
			}
		}
	}

	return true;

}

/*
void HDAG::printHDAG_MLCM(HDAG &G, string filename) {

	ofstream f(filename.c_str());

	f << "\\documentclass{standalone}" << '\n';
	f << "\\usepackage{tikz}" << '\n';
	f << "\\usepackage{amsmath}" << '\n';
	f << "\\usepackage{amsfonts}" << '\n';
	f << "\\usepackage{comment}" << '\n';
	f << "\\begin{document}" << '\n';
	f << "\\begin{tikzpicture}[thick]" << '\n';

	// Print the nodes
	unsigned index_l = 0;

	unsigned max_dim_level = 0;
	double x_position_label_cost = ((G.getLevNumber() - 1) * 15) / 2;

	f << "% VERTICES" << '\n';

	for (unsigned i = 0; i < G.getLevNumber(); ++i) {

		unsigned index_i = 0;

		f << "% Node of level " << i << '\n';

		if (G.getLevel(i).size() > max_dim_level)
			max_dim_level = G.getLevel(i).size();

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			if (G.Os[i][j] == 1) {

				f << "\\node [draw, circle] (" << G.IDs[i][j] << "-" << i
						<< ") at (" << index_l << ", " << index_i << ") {"
						<< G.IDs[i][j] << "};" << '\n';

			} else {

				f << "\\node [draw, circle, dashed] (" << G.IDs[i][j] << "-"
						<< i << ") at (" << index_l << ", " << index_i << ") {"
						<< G.IDs[i][j] << "};" << '\n';

			}

			index_i++;

		}

		f << '\n';

		index_l = index_l + 15;

	}

	f << "% EDGES" << '\n';

	//Print edges
	for (unsigned i = 0; i < G.getLevNumber() - 1; ++i) {

		f << "% Connection between levels " << i << " and " << i + 1 << '\n';

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			for (unsigned k = 0; k < G.getLevel(i)[j].size(); ++k) {

				// Find link (u, v)
				unsigned u = G.IDs[i][j];
				unsigned v = G.LEVELS[i][j][k];

				if (isOriginalLink(u, v, i)) {

					f << "\\draw (" << G.IDs[i][j] << "-" << i << ") -- ("
							<< G.getLevel(i)[j][k] << "-" << i + 1 << ");"
							<< '\n';

				} else {

					f << "\\draw[dashed] (" << G.IDs[i][j] << "-" << i
							<< ") -- (" << G.getLevel(i)[j][k] << "-" << i + 1
							<< ");" << '\n';

				}

			}

		}

	}

	SolutionMLCM S;
	f << "\\node [font=\\Large] at (" << x_position_label_cost << ", " << max_dim_level + 1
			<< ") {MLMC-Solution cost = " << S.getCost(G) << "};" << "'\n";

	f << "\\end{tikzpicture}" << '\n';
	f << "\\end{document}" << '\n';

	f.close();

	return;

}

void HDAG::printHDAG_IMLCM(HDAG &G, string filename) {

	ofstream f(filename.c_str());

	f << "\\documentclass{standalone}" << '\n';
	f << "\\usepackage{tikz}" << '\n';
	f << "\\usepackage{amsmath}" << '\n';
	f << "\\usepackage{amsfonts}" << '\n';
	f << "\\usepackage{comment}" << '\n';
	f << "\\begin{document}" << '\n';
	f << "\\begin{tikzpicture}[thick]" << '\n';

	// Print the nodes
	unsigned index_l = 0;

	unsigned max_dim_level = 0;
	double x_position_label_cost = ((G.getLevNumber() - 1) * 15) / 2;

	f << "% VERTICES" << '\n';

	for (unsigned i = 0; i < G.getLevNumber(); ++i) {

		unsigned index_i = 0;

		f << "% Node of level " << i << '\n';

		if (G.getLevel(i).size() > max_dim_level)
			max_dim_level = G.getLevel(i).size();

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			if (G.Os[i][j] == 1) {

				f << "\\node [draw, circle] (" << G.IDs[i][j] << "-" << i
						<< ") at (" << index_l << ", " << index_i << ") {"
						<< G.IDs[i][j] << "};" << '\n';

			} else {

				if (G.IDs[i][j] == LEVELS[i].size()) {

					f << "\\node [draw, circle, dashed] (" << G.IDs[i][j] << "-"
							<< i << ") at (" << index_l << ", " << index_i << ") { x };" << '\n';

				}
				else {
					f << "\\node [draw, circle, dashed] (" << G.IDs[i][j] << "-"
						<< i << ") at (" << index_l << ", " << index_i << ") {"
						<< G.IDs[i][j] << "};" << '\n';
				}
			}

			index_i++;

		}

		f << '\n';

		index_l = index_l + 15;

	}

	f << "% EDGES" << '\n';

	//Print edges
	for (unsigned i = 0; i < G.getLevNumber() - 1; ++i) {

		f << "% Connection between levels " << i << " and " << i + 1 << '\n';

		for (unsigned j = 0; j < G.getLevel(i).size(); ++j) {

			for (unsigned k = 0; k < G.getLevel(i)[j].size(); ++k) {

				// Find link (u, v)
				unsigned u = G.IDs[i][j];
				unsigned v = G.LEVELS[i][j][k];

				if (isOriginalLink(u, v, i)) {

					f << "\\draw (" << G.IDs[i][j] << "-" << i << ") -- ("
							<< G.getLevel(i)[j][k] << "-" << i + 1 << ");"
							<< '\n';

				} else {

					f << "\\draw[dashed] (" << G.IDs[i][j] << "-" << i
							<< ") -- (" << G.getLevel(i)[j][k] << "-" << i + 1
							<< ");" << '\n';

				}

			}

		}

	}

	SolutionIMLCM S;
	f << "\\node [font=\\Large] at (" << x_position_label_cost << ", " << max_dim_level + 1
			<< ") {CIGDP-Solution cost = " << S.getCost(G) << "};"
			<< "'\n";

	f << "\\end{tikzpicture}" << '\n';
	f << "\\end{document}" << '\n';

	f.close();

	return;

}
*/
