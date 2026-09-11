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

class GLGRASP : public GRASP {
public:
	GLGRASP(HDAG &I, char *argv[]);
	virtual ~GLGRASP();

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

	// Numero maximo de iteracoes sem melhora (eta_max do artigo; 20 nos
	// experimentos de Charytitsch e Nascimento, 2026).
	void set_eta_max(unsigned x) { eta_max = x; }
	unsigned get_eta_max() const { return eta_max; }

private:

	// Valor usado para "sem vizinho posicionado ainda" (G(u) indefinido).
	static const double GL_INF;

	EmbeddingDistances D;

	unsigned eta_max;

	// G[l][i] : menor distancia de embedding entre o vertice cuja posicao
	// ORIGINAL e i no nivel l e os vizinhos dele ja presentes na solucao
	// parcial (Eq. 11). Vale GL_INF enquanto nenhum vizinho tiver sido
	// posicionado.
	vector<vector<double> > G;

	double compute_G(HDAG &I, HDAG &S, unsigned l, unsigned i);
	void init_G(HDAG &I, HDAG &S);
	void update_G_after_insertion(HDAG &I, HDAG &S, unsigned l, unsigned i);

	// Monta a RCL = { v em CL : G(v) <= xi }, com
	// xi = min G + phi (max G - min G)   (Eq. 12).
	void rebuild_RCL(HDAG &I, HDAG &S);

	// Posicao viavel mais proxima de "target" no nivel l.
	unsigned closest_feasible_position(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned target);

	// Step 5 (a)-(c): monta o conjunto P de posicoes candidatas e sorteia uma.
	unsigned choose_position(HDAG &I, HDAG &S, unsigned l, unsigned i, MTRand &rng);

	// Insere o vertice na posicao escolhida, deslocando o bloco para cima
	// quando a posicao ja estiver ocupada.
	void insert_node(HDAG &I, HDAG &S, unsigned l, unsigned i, unsigned new_pos);

};

#endif /* GLGRASP_H_ */
