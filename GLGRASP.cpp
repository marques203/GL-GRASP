/*
 * GLGRASP.cpp
 *
 * Ver GLGRASP.h. As referencias a "Step N", "Eq. (11)" e "Eq. (12)" sao da
 * Secao 4.3 de Charytitsch e Nascimento (2026).
 */

#include "GLGRASP.h"

const double GLGRASP::GL_INF = 1e18;

GLGRASP::GLGRASP(HDAG &I, char *argv[]) {

	set_eta_max(20);   // default do artigo

	position_rule = POS_RULE_FIGURE;
	unanchored_rule = UNANCHORED_MAXDIS;
	feasible_search = SEARCH_REFERENCE;

	setk(atoi(argv[3]));
	setseed(1);
	setalpha(atof(argv[4]));
	setmaxIter(atoi(argv[6]));
	setlocalsearch(argv[5]);
	setTotalInRCL(0);
	setTotalToAssign(0);
	set_best_cost(INT64_MAX);
	set_start_algorithm(chrono::steady_clock::now());
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

	MTRand rng(get_rng_seed());
	set_best_out_file(argv[8]);
	set_complete_out_file(argv[9]);
	setTimeLimit(atof(argv[10]));

	allocate_best_solution(I);

	// Contador de iteracoes sem melhora (it_max no fluxograma do artigo).
	unsigned it_without_improvement = 0;

	for (int i = 0; i < (int) getmaxIter(); i++) {

		// A semente da construcao tem de depender da repeticao, e nao so do
		// indice da iteracao: com "i % 100" a iteracao i usava a mesma semente
		// em todas as repeticoes, tornando-as dependentes entre si (desvios
		// subestimados e p-valores otimistas no Wilcoxon). O codigo de
		// referencia nao semeia nada, logo cada iteracao la e independente.
		//
		// O multiplicador e maior que o max_it maximo aceito pela validacao
		// de entrada (1000000), logo (repeticao, iteracao) -> semente e
		// injetivo: duas repeticoes nunca reaproveitam a mesma semente.
		setseed(get_rng_seed() * 1000003u + (unsigned) i);

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

		long long elapsed = chrono::duration_cast<chrono::microseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000000);

		if (current >= getTimeLimit()) {
			break;
		}

		// Criterio de parada por estagnacao (eta_max do artigo).
		if (it_without_improvement >= get_eta_max()) {
			break;
		}

	}

	long long elapsed = chrono::duration_cast<chrono::microseconds>(
			get_end_iteration() - get_start_algorithm()).count();
	double current = (double) ((double) elapsed / (double) 1000000);

	print_summary_complete_out_file(get_best_cost(), current);
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

		set_end_iteration(chrono::steady_clock::now());

		long long elapsed = chrono::duration_cast<chrono::microseconds>(
				get_end_iteration() - get_start_algorithm()).count();
		double current = (double) ((double) elapsed / (double) 1000000);

		if (best_cost < get_best_cost()) {
			set_best_cost(best_cost);
			set_time_to_best(current);
		}

		print_complete_data_out_file(best_cost, best_cost, current);

		delete_all_iterative_structure(I);

	}

	return;

}

void GLGRASP::construction(HDAG &I, HDAG &S) {

	initialize_all_iteration_structure(I);

	set_start_iteration(chrono::steady_clock::now());

	MTRand rng(getseed());

	// Step 1 e 2: I' recebe apenas os vertices originais; CL fica com os
	// incrementais.
	allocate_solution(I, S);

	// Step 3: G(u) para todo u em CL (Eq. 11).
	init_G(I, S);

	// Step 4: RCL a partir de xi (Eq. 12).
	rebuild_RCL(I, S);

	while (getTotalToAssign() > 0) {

		// Step 5: escolhe aleatoriamente um vertice da RCL, com probabilidade
		// uniforme sobre a RCL inteira.
		pair<unsigned, unsigned> p = extract_from_RCL_uniform(rng);

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

double GLGRASP::effective_G(unsigned l, unsigned i) const {

	double g = G[l][i];

	// Semantica do codigo de referencia: o vertice sem vizinho posicionado
	// recebe a maior distancia do grafo, valor que PARTICIPA do max G da
	// Eq. (12) e portanto infla xi.
	if (g >= GL_INF && unanchored_rule == UNANCHORED_MAXDIS)
		return D.max_distance();

	return g;

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
			if (CL[l][i] == 1) {
				double g = effective_G(l, i);
				if (g < GL_INF) {
					has_finite = true;
					if (g < g_min) g_min = g;
					if (g > g_max) g_max = g;
				}
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

	// Com UNANCHORED_INF, vertices com G indefinido ficam de fora enquanto
	// existir candidato com distancia finita: eles sao os ultimos a serem
	// posicionados. Com UNANCHORED_MAXDIS eles valem max_dis = g_max e so
	// entram se phi == 1, como no codigo de referencia.
	for (unsigned l = 0; l < G.size(); ++l)
		for (unsigned i = 0; i < G[l].size(); ++i)
			if (CL[l][i] == 1) {
				double g = effective_G(l, i);
				if (g < GL_INF && g <= xi)
					add_in_RCL(l, i);
			}

	// Rede de seguranca: xi >= g_min garante ao menos um candidato, mas uma
	// eventual imprecisao de ponto flutuante nao pode deixar a RCL vazia
	// (extract_from_RCL entraria em laco infinito).
	if (getTotalInRCL() == 0) {
		for (unsigned l = 0; l < G.size(); ++l)
			for (unsigned i = 0; i < G[l].size(); ++i)
				if (CL[l][i] == 1 && effective_G(l, i) == g_min) {
					add_in_RCL(l, i);
					return;
				}
	}

	return;

}

pair<unsigned, unsigned> gl_rcl_pick_uniform(
		vector<vector<pair<unsigned, unsigned> > > &RCL,
		unsigned total, MTRand &rng) {

	// rebuild_RCL garante RCL nao vazia; a guarda existe so para nao haver
	// leitura fora de faixa se essa invariante for quebrada no futuro.
	if (total == 0)
		return make_pair(0u, 0u);

	unsigned r = rng.randInt(total - 1);

	for (unsigned l = 0; l < RCL.size(); ++l) {

		if (r < RCL[l].size()) {

			pair<unsigned, unsigned> tmp = RCL[l][r];
			RCL[l].erase(RCL[l].begin() + r);

			return tmp;

		}

		r -= RCL[l].size();

	}

	// Inalcancavel enquanto total == soma dos tamanhos.
	return make_pair(0u, 0u);

}

pair<unsigned, unsigned> GLGRASP::extract_from_RCL_uniform(MTRand &rng) {

	pair<unsigned, unsigned> tmp =
			gl_rcl_pick_uniform(getRCL(), getTotalInRCL(), rng);

	if (getTotalInRCL() > 0)
		getTotalInRCL()--;

	return tmp;

}

unsigned GLGRASP::feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i,
		unsigned target) {

	if (feasible_search == SEARCH_REFERENCE)
		return reference_feasible_position(I, S, l, i, target);

	return closest_feasible_position(I, S, l, i, target);

}

bool GLGRASP::can_shift_up(HDAG &I, HDAG &S, unsigned l, unsigned pos) {

	auto k = getk();

	auto &Pos = I.getPos();
	auto &IDs_S = S.getIDs();
	auto &Os_S = S.getOs();

	if (pos == 0)
		return false;

	// Vertice incremental nao tem restricao de deslocamento.
	if (!Os_S[l][pos - 1])
		return true;

	// Original: Claim 2, |pi_0(v) - pi(v)| <= k, nos dois limites.
	unsigned id_ = IDs_S[l][pos - 1];
	unsigned origin_pos = Pos[l][id_];

	if (pos > origin_pos + k)
		return false;

	if (origin_pos > k && pos < origin_pos - k)
		return false;

	return true;

}

unsigned GLGRASP::reference_feasible_position(HDAG &I, HDAG &S, unsigned l,
		unsigned i, unsigned target) {

	auto &LEVELS_S = S.getLEVELS();

	unsigned n = LEVELS_S[l].size();

	// ind do g_vq: primeiro slot livre do nivel.
	unsigned first_free = n;
	for (unsigned p = 0; p < n; ++p) {
		if (is_position_free(S, l, p)) {
			first_free = p;
			break;
		}
	}

	// A construcao so chama isto enquanto ha vertice a posicionar, logo
	// sempre existe slot livre; a guarda evita indice invalido se a
	// invariante for quebrada no futuro.
	if (first_free == n)
		return target < n ? target : n - 1;

	// "if q3 >= ind: pos_final = ind" -- alvo no slot livre ou alem dele.
	if (target >= first_free)
		return first_free;

	// Caminha do slot livre PARA BAIXO em direcao ao alvo, perguntando a cada
	// passo se o vertice que esta em p-1 pode subir para p. Para no primeiro
	// passo bloqueado.
	unsigned con = 0;
	for (unsigned p = first_free; p > target; --p) {
		if (!can_shift_up(I, S, l, p))
			break;
		++con;
	}

	// "if con == ind - q3: pos_final = q3  else: pos_final = ind - con"
	if (con == first_free - target)
		return target;

	return first_free - con;

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

unsigned GLGRASP::position_without_anchor(HDAG &I, HDAG &S, unsigned l,
		unsigned i, MTRand &rng) {

	// Codigo de referencia: random.randint(0, n - 1).
	if (unanchored_rule == UNANCHORED_MAXDIS) {
		unsigned n = S.getLEVELS()[l].size();
		unsigned target = rng.randInt(n - 1);
		return feasible_position(I, S, l, i, target);
	}

	// Nossa semantica: mantem o vertice o mais perto possivel da posicao que
	// ele ocupa na instancia.
	return feasible_position(I, S, l, i, i);

}

unsigned GLGRASP::choose_position(HDAG &I, HDAG &S, unsigned l, unsigned i,
		MTRand &rng) {

	auto &IDs = I.getIDs();
	auto &Pos = I.getPos();
	auto &Pos_S = S.getPos();
	auto &CL = getCL();

	unsigned id_u = IDs[l][i];

	// Alvo escolhido em cada nivel vizinho (l+1 e l-1), antes de qualquer
	// ajuste de viabilidade.
	vector<unsigned> side_target;

	// Candidatas ja ajustadas para viabilidade (so a leitura POS_RULE_TEXT).
	vector<unsigned> P;

	// Step 5 (a): em cada nivel vizinho que tenha algum vizinho de u ja
	// posicionado, localiza o vizinho MAIS PROXIMO e o MAIS DISTANTE (em
	// distancia de embedding) e define o alvo daquele lado.
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

		unsigned pos_min = Pos_S[lv][id_min];
		unsigned pos_max = Pos_S[lv][id_max];

		if (position_rule == POS_RULE_FIGURE) {

			// Figura 3 / func_ale: uniforme entre a posicao do vizinho mais
			// proximo, a do mais distante e o ponto medio das duas. O ponto
			// medio DENTRO de um mesmo nivel so existe nesta leitura.
			side_target.push_back(gl_func_ale(pos_min, pos_max, rng));

		} else {

			// Texto do Step 5: sorteia entre mais proximo e mais distante, e
			// ja converte para a posicao viavel mais proxima no nivel l.
			unsigned chosen = (rng.randInt(1) == 0) ? pos_min : pos_max;
			P.push_back(feasible_position(I, S, l, i, chosen));

		}

	}

	if (position_rule == POS_RULE_FIGURE) {

		if (side_target.size() == 0)
			return position_without_anchor(I, S, l, i, rng);

		// Entre niveis, func_ale de novo; com um unico nivel vizinho o alvo
		// daquele lado e usado diretamente.
		unsigned target = (side_target.size() == 2)
				? gl_func_ale(side_target[0], side_target[1], rng)
				: side_target[0];

		// O ajuste de viabilidade e aplicado UMA vez, sobre o alvo final.
		return feasible_position(I, S, l, i, target);

	}

	// Step 5 (b): acrescenta a posicao viavel mais proxima da media das que
	// ja estao em P.
	if (P.size() > 0) {
		unsigned long long sum = 0;
		for (unsigned j = 0; j < P.size(); ++j)
			sum += P[j];
		unsigned avg = (unsigned) (sum / P.size());
		P.push_back(feasible_position(I, S, l, i, avg));
	}

	if (P.size() == 0)
		return position_without_anchor(I, S, l, i, rng);

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

			// Claim 2 e |pi_0(v) - pi(v)| <= k, ou seja, tem DOIS limites. O
			// deslocamento aqui e sempre para cima, logo este teste nunca
			// deveria falhar; e a mesma verificacao defensiva que o
			// check_position do codigo de referencia faz.
			if (origin_pos > k && potential_new_pos < origin_pos - k)
				return 1;

		}

		last_in_block++;

	}

	// O bloco inteiro nao pode ser deslocado
	if (last_in_block == n_nodes_in_lev)
		return 1;

	return 0;

}
