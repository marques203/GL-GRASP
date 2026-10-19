/*
 * GLGRASP.h
 *
 * Graph Learning GRASP (GL-GRASP) para o C-IGDP, seguindo a metodologia de
 * Charytitsch e Nascimento (2026), Secao 4.3.
 *
 * A diferenca em relacao ao GRASP_v3 (heuristica C3 de Napoletano et al.,
 * 2019) esta apenas na FASE DE CONSTRUCAO:
 *
 *   - o criterio guloso deixa de ser "numero de cruzamentos adicionais" e
 *     passa a ser a distancia estrutural G(u) obtida dos node embeddings
 *     (Eq. 11), lida de um arquivo gerado pelo modulo Python;
 *   - a posicao de insercao deixa de ser a que minimiza cruzamentos e passa
 *     a ser sorteada entre as candidatas do Step 5 (vizinho mais proximo,
 *     mais distante e media).
 *
 * A busca local (swap + insert) e reaproveitada sem alteracao da classe base,
 * como o artigo especifica.
 */

#include "GRASP.h"
#include "EmbeddingDistances.h"

#ifndef GLGRASP_H_
#define GLGRASP_H_

using namespace std;

// As duas funcoes abaixo concentram a logica nova do Step 5 e do sorteio na
// RCL. Sao livres (e nao metodos privados) para que os testes de unidade em
// tests/test_glgrasp_units.cpp exercitem exatamente o mesmo codigo que roda
// em producao, sem precisar instanciar um GLGRASP inteiro.

// func_ale do Cb_*.py: uniforme entre a, b e o ponto medio de a e b.
inline unsigned gl_func_ale(unsigned a, unsigned b, MTRand &rng) {
	// unsigned long tem 32 bits no MinGW, logo a soma precisa ser feita em
	// 64 bits para o ponto medio nao transbordar.
	unsigned opts[3] = { a, b, (unsigned) (((unsigned long long) a + b) / 2) };
	return opts[rng.randInt(2)];
}

// Sorteia uniformemente UMA entrada da RCL achatada e a remove. "total" tem
// de ser a soma dos tamanhos dos niveis. Devolve (nivel, indice) do vertice.
// Contraste com GRASP::extract_from_RCL, que sorteia o nivel primeiro: ali
// P(v) = 1/(L * |RCL_l|), aqui P(v) = 1/|RCL|.
pair<unsigned, unsigned> gl_rcl_pick_uniform(
		vector<vector<pair<unsigned, unsigned> > > &RCL,
		unsigned total, MTRand &rng);

class GLGRASP : public GRASP {
public:

	// O artigo e internamente inconsistente sobre o Step 5: o TEXTO descreve
	// a escolha entre vizinho mais proximo e mais distante por nivel, seguida
	// da media do conjunto P; a FIGURA 3 (e o codigo de referencia Cb_*.py)
	// sorteia uniformemente entre tres alvos, {mais proximo, mais distante,
	// ponto medio}, por nivel, e de novo entre os niveis. As duas leituras
	// estao implementadas.
	enum PositionRule {
		POS_RULE_FIGURE,   // func_ale do Cb_*.py (default)
		POS_RULE_TEXT      // leitura literal do Step 5
	};

	// Tratamento de G(u) indefinido (nenhum vizinho de u posicionado ainda).
	// INF: u fica fora da RCL enquanto houver candidato com distancia finita,
	//      e o intervalo [min G, max G] da Eq. (12) ignora esses vertices.
	// MAXDIS: G(u) = maior distancia do grafo, como no Cb_*.py; o valor entra
	//      no max G da Eq. (12) e a posicao de u e sorteada uniformemente.
	enum UnanchoredRule {
		UNANCHORED_MAXDIS, // reproduz o codigo de referencia (default)
		UNANCHORED_INF
	};

	// Como achar uma posicao viavel quando o alvo do Step 5 nao serve.
	// REFERENCE: g_vq do Cb_*.py. Se o alvo esta no slot livre ou alem dele,
	//   insere no primeiro slot livre; senao caminha do slot livre PARA BAIXO
	//   em direcao ao alvo, parando onde a Claim 2 bloqueia. A posicao final
	//   fica sempre entre o alvo e o primeiro slot livre, nunca abaixo do alvo.
	// SPIRAL: busca em espiral em torno do alvo, aceitando posicao dos dois
	//   lados; e a leitura literal de "posicao viavel mais proxima".
	enum FeasibleSearch {
		SEARCH_REFERENCE,  // (default)
		SEARCH_SPIRAL
	};

	GLGRASP(HDAG &I, char *argv[]);
	virtual ~GLGRASP();

	void set_position_rule(PositionRule r) { position_rule = r; };
	PositionRule get_position_rule() const { return position_rule; };

	void set_unanchored_rule(UnanchoredRule r) { unanchored_rule = r; };
	UnanchoredRule get_unanchored_rule() const { return unanchored_rule; };

	void set_feasible_search(FeasibleSearch s) { feasible_search = s; };
	FeasibleSearch get_feasible_search() const { return feasible_search; };

	virtual void algorithm(HDAG &I, int argc, char *argv[]);
	virtual void run(HDAG &I, HDAG &S);
	virtual void construction(HDAG &I, HDAG &S);
	virtual void allocate_solution(HDAG &I, HDAG &S);
	virtual unsigned check_violation(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned pos);

	// Nao utilizados nesta heuristica (a insercao e feita por insert_node),
	// mas exigidos pela interface da classe base.
	virtual void add_node_to_solution(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng) {};
	virtual unsigned find_feasible_position_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_original_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };
	virtual unsigned find_feasible_position_incremental_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned bc) { return 0; };

	// Carrega o arquivo de distancias produzido pelo modulo Python.
	bool load_distances(const string &filename, HDAG &I) { return D.load(filename, I); }

private:

	// Valor usado para "sem vizinho posicionado ainda" (G(u) indefinido).
	static const double GL_INF;

	PositionRule position_rule;
	UnanchoredRule unanchored_rule;
	FeasibleSearch feasible_search;

	EmbeddingDistances D;

	// G[l][i] : menor distancia de embedding entre o vertice cuja posicao
	// ORIGINAL e i no nivel l e os vizinhos dele ja presentes na solucao
	// parcial (Eq. 11). Vale GL_INF enquanto nenhum vizinho tiver sido
	// posicionado.
	vector<vector<double> > G;

	double compute_G(HDAG &I, HDAG &S, unsigned l, unsigned i);
	void init_G(HDAG &I, HDAG &S);
	void update_G_after_insertion(HDAG &I, HDAG &S, unsigned l, unsigned i);

	// G(l, i) como entra na Eq. (12): resolve G indefinido conforme a regra
	// escolhida. Devolve GL_INF quando o vertice deve ser ignorado.
	double effective_G(unsigned l, unsigned i) const;

	// Monta a RCL = { v em CL : G(v) <= xi }, com
	// xi = min G + phi (max G - min G)   (Eq. 12).
	void rebuild_RCL(HDAG &I, HDAG &S);

	// Sorteia UNIFORMEMENTE entre todos os vertices da RCL, como o Step 5 pede
	// ("select a vertex v in RCL at random") e como faz o codigo de referencia.
	// Difere do GRASP::extract_from_RCL da classe base, que sorteia primeiro o
	// nivel e depois o indice dentro dele -- isso superamostra niveis com
	// poucos candidatos. A base fica intacta para GRASP1/2/3.
	pair<unsigned, unsigned> extract_from_RCL_uniform(MTRand &rng);

	// Posicao viavel para o vertice (l, i) dado o alvo do Step 5; despacha
	// para a regra escolhida em feasible_search.
	unsigned feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned target);

	// SEARCH_SPIRAL: posicao viavel mais proxima de "target", dos dois lados.
	unsigned closest_feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned target);

	// SEARCH_REFERENCE: g_vq do Cb_*.py.
	unsigned reference_feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned target);

	// O vertice que hoje ocupa a posicao pos-1 pode subir para pos sem violar
	// a Claim 2? (check_position do codigo de referencia.)
	bool can_shift_up(HDAG &I, HDAG &S, unsigned l, unsigned pos);

	// Step 5 (a)-(c): monta as posicoes candidatas e sorteia uma.
	unsigned choose_position(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);

	// Posicao do vertice sem nenhum vizinho posicionado.
	unsigned position_without_anchor(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);

	// Insere o vertice na posicao escolhida, deslocando o bloco para cima
	// quando a posicao ja estiver ocupada.
	void insert_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned new_pos);

};

#endif /* GLGRASP_H_ */
