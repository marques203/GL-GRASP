/*
 * EmbeddingDistances.h
 *
 * Distancias estruturais por arco, geradas pelo modulo Python
 * (embeddings/build_distances.py) a partir dos node embeddings.
 *
 * Formato do arquivo lido:
 *     n_arcos
 *     nivel_u id_u nivel_v id_v distancia
 *     ...
 * com nivel_v == nivel_u + 1 e os ids LOCAIS ao nivel, exatamente como no
 * arquivo de instancia lido por HDAG::read_instance.
 */

#ifndef EMBEDDINGDISTANCES_H_
#define EMBEDDINGDISTANCES_H_

#include <string>
#include <vector>
#include <utility>
#include "HDAG.h"

class EmbeddingDistances {
public:
	EmbeddingDistances() : n_arcs(0) {}

	// Le o arquivo e indexa as distancias nos dois sentidos.
	// Devolve false (e explica no cerr) se o arquivo nao existir ou nao
	// for compativel com a instancia I.
	bool load(const std::string &filename, HDAG &I);

	// Vizinhos de (l, id) no nivel l+1: pares (id_do_vizinho, distancia).
	const std::vector<std::pair<unsigned, double> > &forward(unsigned l, unsigned id) const {
		return FWD[l][id];
	}

	// Vizinhos de (l, id) no nivel l-1: pares (id_do_vizinho, distancia).
	const std::vector<std::pair<unsigned, double> > &backward(unsigned l, unsigned id) const {
		return BWD[l][id];
	}

	unsigned get_n_arcs() const { return n_arcs; }

private:
	unsigned n_arcs;

	// FWD[l][id] : arcos que saem de (l, id) para o nivel l+1
	std::vector<std::vector<std::vector<std::pair<unsigned, double> > > > FWD;

	// BWD[l][id] : arcos que chegam em (l, id) vindos do nivel l-1
	std::vector<std::vector<std::vector<std::pair<unsigned, double> > > > BWD;

};

#endif /* EMBEDDINGDISTANCES_H_ */
