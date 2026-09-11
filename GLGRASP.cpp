/*
 * GLGRASP.cpp
 *
 * Ver GLGRASP.h. As referencias a "Step N", "Eq. (11)" e "Eq. (12)" sao da
 * Secao 4.3 de Charytitsch e Nascimento (2026).
 */

#include "GLGRASP.h"

const double GLGRASP::GL_INF = 1e18;

GLGRASP::GLGRASP(HDAG &I, char *argv[]) : eta_max(20) {

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

GLGRASP::~GLGRASP() {}

void GLGRASP::algorithm(HDAG &I, int argc, char *argv[]) {

	SolutionIMLCM s;
	set_best_cost(s.getCost(I));
	set_time_to_best(0.001);
	cout << "instance cost " << get_best_cost();

	double initial_alpha = getalpha();

	MTRand rng;
	set_best_out_file(argv[8]);
	set_complete_out_file(argv[9]);
	setTimeLimit(atof(argv[10]));

	allocate_best_solution(I);

	// Contador de iteracoes sem melhora (it_max no fluxograma do artigo).
	unsigned it_without_improvement = 0;

	for (int i = 0; i < (int) getmaxIter(); i++) {

		setseed(i % 100);

		// phi e sorteado a cada iteracao do GRASP (Eq. 12); alpha = 0 na
		// linha de comando ativa esse modo, como no GRASP_v3.
		if (initial_alpha == 0) {
			setalpha(rng.randDblExc(1));
		}

		double cost_before = get_best_cost();

		// Mesmo padrao de alocacao do GRASP_v3, para que as duas heuristicas
		// tenham o mesmo comportamento de memoria na comparacao.
		HDAG *S = new HDAG;
		run(I, *S);
		S = nullptr;

		if (get_best_cost() < cost_before) {
			it_without_improvement = 0;
		} else {
			it_without_improvement++;
		}

		int elapsed = chrono::duration_cast<chrono::milliseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000);

		if (current >= getTimeLimit()) {
			break;
		}

		// Criterio de parada por estagnacao (eta_max do artigo).
		if (it_without_improvement >= eta_max) {
			break;
		}

	}

	int elapsed = chrono::duration_cast<chrono::milliseconds>(
			get_end_iteration() - get_start_algorithm()).count();
	double current = (double) ((double) elapsed / (double) 1000);

	print_complete_data_out_file(get_best_cost(), current);
	print_best_data_out_file();
	close_out_file();

	return;

}

void GLGRASP::run(HDAG &I, HDAG &S) {

	this->construction(I, S);

	if (strcmp(getlocalsearch(), "no") != 0) {

		// Busca local identica a de Napoletano et al. (2019), como o artigo
		// especifica: nao ha nada de embedding aqui.
		this->localsearch(I, S);

	} else {

		SolutionIMLCM s;
		s.initialize_costsForLevels(I.getLEVELS().size());

		double best_cost = s.getCost(S);

		S.set_cost(best_cost);

		set_end_iteration(chrono::system_clock::now());

		int elapsed = chrono::duration_cast<chrono::milliseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000);

		if (best_cost < get_best_cost()) {
			set_best_cost(best_cost);
			set_time_to_best(current);
		}

		print_complete_data_out_file(best_cost, current);

		delete_all_iterative_structure(I);

	}

	return;

}

void GLGRASP::construction(HDAG &I, HDAG &S) {

	initialize_all_iteration_structure(I);

	set_start_iteration(chrono::system_clock::now());

	MTRand rng(getseed());

	// Step 1 e 2: I' recebe apenas os vertices originais; CL fica com os
	// incrementais.
	allocate_solution(I, S);

	// Step 3: G(u) para todo u em CL (Eq. 11).
	init_G(I, S);

	// Step 4: RCL a partir de xi (Eq. 12).
	rebuild_RCL(I, S);

	while (getTotalToAssign() > 0) {

		// Step 5: escolhe aleatoriamente um vertice da RCL.
		pair<unsigned, unsigned> p = extract_from_RCL(rng);

		// Step 5 (a)-(c): monta as posicoes candidatas e sorteia uma.
		unsigned new_pos = choose_position(I, S, p.first, p.second, rng);

		insert_node(I, S, p.first, p.second, new_pos);

		// Step 6 (com a otimizacao descrita no artigo): so os vizinhos do
		// vertice recem-inserido precisam ter G reavaliado.
		update_G_after_insertion(I, S, p.first, p.second);

		rebuild_RCL(I, S);

	}

	if (!S.is_feasible(I, getk()))
		cerr << " NOT FEASIBLE!!! " << endl;

	return;

}

double GLGRASP::compute_G(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &CL = getCL();

	unsigned id_u = IDs[l][i];

	double best = GL_INF;

	// Vizinhos no nivel seguinte.
	if (l + 1 < G.size()) {
		const vector<pair<unsigned, double> > &nbs = D.forward(l, id_u);
		for (unsigned j = 0; j < nbs.size(); ++j) {
			unsigned pos_v = Pos[l + 1][nbs[j].first];
			// CL == 0 significa "ja esta na solucao parcial".
			if (CL[l + 1][pos_v] == 0 && nbs[j].second < best)
				best = nbs[j].second;
		}
	}

	// Vizinhos no nivel anterior.
	if (l > 0) {
		const vector<pair<unsigned, double> > &nbs = D.backward(l, id_u);
		for (unsigned j = 0; j < nbs.size(); ++j) {
			unsigned pos_v = Pos[l - 1][nbs[j].first];
			if (CL[l - 1][pos_v] == 0 && nbs[j].second < best)
				best = nbs[j].second;
		}
	}

	return best;

}

void GLGRASP::init_G(HDAG &I, HDAG &S) {

	auto &LEVELS = I.getLEVELS();
	auto &CL = getCL();

	G.clear();
	G.resize(LEVELS.size());

	for (unsigned l = 0; l < LEVELS.size(); ++l) {

		G[l].assign(LEVELS[l].size(), GL_INF);

		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			if (CL[l][i] == 1)
				G[l][i] = compute_G(I, S, l, i);
		}

	}

	return;

}

void GLGRASP::update_G_after_insertion(HDAG &I, HDAG &S, unsigned l, unsigned i) {

	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &CL = getCL();

	unsigned id_u = IDs[l][i];

	// O vertice (l, i) acabou de entrar na solucao: para cada vizinho ainda
	// nao posicionado, G so pode diminuir, e no maximo ate a distancia deste
	// novo arco. Nao ha necessidade de recalcular o minimo do zero.
	if (l + 1 < G.size()) {
		const vector<pair<unsigned, double> > &nbs = D.forward(l, id_u);
		for (unsigned j = 0; j < nbs.size(); ++j) {
			unsigned pos_v = Pos[l + 1][nbs[j].first];
			if (CL[l + 1][pos_v] == 1 && nbs[j].second < G[l + 1][pos_v])
				G[l + 1][pos_v] = nbs[j].second;
		}
	}

	if (l > 0) {
		const vector<pair<unsigned, double> > &nbs = D.backward(l, id_u);
		for (unsigned j = 0; j < nbs.size(); ++j) {
			unsigned pos_v = Pos[l - 1][nbs[j].first];
			if (CL[l - 1][pos_v] == 1 && nbs[j].second < G[l - 1][pos_v])
				G[l - 1][pos_v] = nbs[j].second;
		}
	}

	return;

}

void GLGRASP::rebuild_RCL(HDAG &I, HDAG &S) {

	auto &RCL = getRCL();
	auto &CL = getCL();

	for (unsigned l = 0; l < RCL.size(); ++l)
		RCL[l].clear();

	setTotalInRCL(0);

	double g_min = GL_INF;
	double g_max = -1.0;
	bool has_finite = false;

	for (unsigned l = 0; l < G.size(); ++l) {
		for (unsigned i = 0; i < G[l].size(); ++i) {
			if (CL[l][i] == 1 && G[l][i] < GL_INF) {
				has_finite = true;
				if (G[l][i] < g_min) g_min = G[l][i];
				if (G[l][i] > g_max) g_max = G[l][i];
			}
		}
	}

	// Nenhum candidato tem vizinho ja posicionado (pode acontecer quando um
	// componente e inteiramente incremental). Nesse caso todos entram na RCL,
	// para que a construcao nao trave.
	if (!has_finite) {
		for (unsigned l = 0; l < G.size(); ++l)
			for (unsigned i = 0; i < G[l].size(); ++i)
				if (CL[l][i] == 1)
					add_in_RCL(l, i);
		return;
	}

	// Eq. (12): xi = min G + phi (max G - min G).
	double xi = g_min + getalpha() * (g_max - g_min);

	// Vertices com G indefinido ficam de fora enquanto existir candidato com
	// distancia finita: eles sao os ultimos a serem posicionados.
	for (unsigned l = 0; l < G.size(); ++l)
		for (unsigned i = 0; i < G[l].size(); ++i)
			if (CL[l][i] == 1 && G[l][i] < GL_INF && G[l][i] <= xi)
				add_in_RCL(l, i);

	// Rede de seguranca: xi >= g_min garante ao menos um candidato, mas uma
	// eventual imprecisao de ponto flutuante nao pode deixar a RCL vazia
	// (extract_from_RCL entraria em laco infinito).
	if (getTotalInRCL() == 0) {
		for (unsigned l = 0; l < G.size(); ++l)
			for (unsigned i = 0; i < G[l].size(); ++i)
				if (CL[l][i] == 1 && G[l][i] == g_min) {
					add_in_RCL(l, i);
					return;
				}
	}

	return;

}

unsigned GLGRASP::closest_feasible_position(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned target) {

	auto &LEVELS_S = S.getLEVELS();

	unsigned n = LEVELS_S[l].size();

	if (target >= n)
		target = n - 1;

	if (check_violation(I, S, l, i, target) == 0)
		return target;

	for (unsigned delta = 1; delta < n; ++delta) {

		if (target >= delta) {
			unsigned p = target - delta;
			if (check_violation(I, S, l, i, p) == 0)
				return p;
		}

		unsigned p = target + delta;
		if (p < n && check_violation(I, S, l, i, p) == 0)
			return p;

	}

	// Ultimo recurso: qualquer posicao livre.
	for (unsigned p = 0; p < n; ++p)
		if (is_position_free(S, l, p))
			return p;

	return target;

}

unsigned GLGRASP::choose_position(HDAG &I, HDAG &S, unsigned l, unsigned i,
		MTRand &rng) {

	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &Pos_S = S.getPos();
	auto &CL = getCL();

	unsigned id_u = IDs[l][i];

	vector<unsigned> P;

	// Step 5 (a): para cada nivel vizinho que tenha algum vizinho de u ja
	// posicionado, sorteia entre a posicao do vizinho MAIS PROXIMO e a do
	// MAIS DISTANTE (em distancia de embedding).
	for (int side = 0; side < 2; ++side) {

		unsigned lv;
		const vector<pair<unsigned, double> > *nbs;

		if (side == 0) {
			if (l + 1 >= G.size()) continue;
			lv = l + 1;
			nbs = &D.forward(l, id_u);
		} else {
			if (l == 0) continue;
			lv = l - 1;
			nbs = &D.backward(l, id_u);
		}

		bool found = false;
		double d_min = 0, d_max = 0;
		unsigned id_min = 0, id_max = 0;

		for (unsigned j = 0; j < nbs->size(); ++j) {

			unsigned id_v = (*nbs)[j].first;
			double d = (*nbs)[j].second;

			// so vizinhos ja presentes na solucao parcial
			if (CL[lv][Pos[lv][id_v]] != 0)
				continue;

			if (!found) {
				found = true;
				d_min = d_max = d;
				id_min = id_max = id_v;
			} else {
				if (d < d_min) { d_min = d; id_min = id_v; }
				if (d > d_max) { d_max = d; id_max = id_v; }
			}

		}

		if (!found)
			continue;

		unsigned chosen_id = (rng.randInt(1) == 0) ? id_min : id_max;

		// A posicao do vizinho escolhido, no nivel vizinho, e o alvo; a
		// candidata e a posicao viavel mais proxima dele no nivel l.
		P.push_back(closest_feasible_position(I, S, l, i, Pos_S[lv][chosen_id]));

	}

	// Step 5 (b): acrescenta a posicao viavel mais proxima da media das que
	// ja estao em P.
	if (P.size() > 0) {
		unsigned long sum = 0;
		for (unsigned j = 0; j < P.size(); ++j)
			sum += P[j];
		unsigned avg = (unsigned) (sum / P.size());
		P.push_back(closest_feasible_position(I, S, l, i, avg));
	}

	// Nenhum vizinho posicionado (G(u) indefinido): mantem o vertice o mais
	// perto possivel da posicao que ele ocupa na instancia.
	if (P.size() == 0)
		return closest_feasible_position(I, S, l, i, i);

	// Step 5 (c): escolhe uma das candidatas ao acaso.
	return P[rng.randInt(P.size() - 1)];

}

void GLGRASP::insert_node(HDAG &I, HDAG &S, unsigned l, unsigned i,
		unsigned new_pos) {

	auto &IDs = I.getIDs();
	auto &Os = I.getOs();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Pos_S = S.getPos();
	auto &Os_S = S.getOs();

	auto &CL = getCL();

	auto &upmost_inserted_node = getupmost_inserted_node();

	unsigned n_nodes_in_level = LEVELS_S[l].size();

	if (!is_position_free(S, l, new_pos)) {

		// Encontra a primeira posicao livre acima de new_pos e desloca o
		// bloco para cima, abrindo espaco.
		unsigned i_ = new_pos;
		while (true) {
			i_++;
			if (is_position_free(S, l, i_))
				break;
		}

		unsigned first_free = i_;
		unsigned first_busy = first_free - 1;

		for (i_ = first_busy; i_ >= new_pos; --i_) {

			unsigned id_i_ = IDs_S[l][i_];
			Pos_S[l][id_i_] = i_ + 1;
			LEVELS_S[l][i_ + 1].swap(LEVELS_S[l][i_]);
			Os_S[l][i_ + 1] = Os_S[l][i_];
			IDs_S[l][i_ + 1] = id_i_;

			// esvazia a posicao liberada
			Os_S[l][i_] = n_nodes_in_level;
			IDs_S[l][i_] = n_nodes_in_level;

			if (is_position_free(S, l, new_pos)) break;

			if (i_ == 0)
				break;

		}

	}

	if (new_pos > upmost_inserted_node[l]) {
		upmost_inserted_node[l] = new_pos;
	}

	unsigned id_u = IDs[l][i];
	IDs_S[l][new_pos] = id_u;
	Pos_S[l][id_u] = new_pos;
	Os_S[l][new_pos] = Os[l][i];
	add_link(I, S, l, i, new_pos);

	// Remove de CL
	CL[l][i] = 0;

	gettoAssign()[l]--;
	getTotalToAssign()--;

	return;

}

void GLGRASP::allocate_solution(HDAG &I, HDAG &S) {

	auto k = getk();
	auto &CL = getCL();
	auto &maxdeg = getmaxdeg();
	auto &toAssign = gettoAssign();
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

				// Remove de CL
				CL[l][i] = 0;

				S.increase_total_nodes(1);

				if (upmost_inserted_node[l] == vn) {
					upmost_inserted_node[l] = i;
				} else if (upmost_inserted_node[l] < i) {
					upmost_inserted_node[l] = i;
				}

			} else {

				maxdeg = max(OD[l][i], maxdeg);

				toAssign[l]++;
				totalToAssign++;

			}

		}

	}

	return;

}

unsigned GLGRASP::check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i,
		unsigned pos) {

	auto k = getk();

	auto &Pos = I.getPos();

	auto &LEVELS_S = S.getLEVELS();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	unsigned n_nodes_in_lev = LEVELS_S[l].size();

	if (is_position_free(S, l, pos))
		return 0;

	// Verifica se o bloco ocupado a partir de "pos" pode ser deslocado uma
	// posicao para cima sem violar a restricao de deslocamento maximo (k)
	// de nenhum vertice original.
	unsigned last_in_block = pos;
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

	// O bloco inteiro nao pode ser deslocado
	if (last_in_block == n_nodes_in_lev)
		return 1;

	return 0;

}
