/*
 * EmbeddingDistances.cpp
 */

#include "EmbeddingDistances.h"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

bool EmbeddingDistances::load(const string &filename, HDAG &I) {

	ifstream f(filename.c_str());

	if (!f.is_open()) {
		cerr << " Nao foi possivel abrir o arquivo de distancias: "
			 << filename << endl;
		return false;
	}

	auto &LEVELS = I.getLEVELS();
	unsigned ln = LEVELS.size();

	FWD.clear();
	BWD.clear();
	FWD.resize(ln);
	BWD.resize(ln);
	for (unsigned l = 0; l < ln; ++l) {
		FWD[l].resize(LEVELS[l].size());
		BWD[l].resize(LEVELS[l].size());
	}

	string line;

	if (!getline(f, line)) {
		cerr << " Arquivo de distancias vazio: " << filename << endl;
		return false;
	}

	unsigned declared = 0;
	{
		istringstream iss(line);
		iss >> declared;
	}

	unsigned read = 0;

	while (getline(f, line)) {

		if (line.find_first_not_of(" \t\r\n") == string::npos)
			continue;

		istringstream iss(line);
		unsigned lu, iu, lv, iv;
		double d;

		if (!(iss >> lu >> iu >> lv >> iv >> d)) {
			cerr << " Linha malformada no arquivo de distancias: " << line << endl;
			return false;
		}

		// A instancia e o arquivo de distancias precisam ser da mesma
		// instancia: sem essa checagem, um arquivo trocado corromperia a
		// construcao silenciosamente.
		if (lu >= ln || lv >= ln) {
			cerr << " Nivel fora da faixa no arquivo de distancias (" << lu
				 << " -> " << lv << "), a instancia tem " << ln << " niveis."
				 << " O arquivo de distancias corresponde a esta instancia?" << endl;
			return false;
		}

		if (lv != lu + 1) {
			cerr << " Arco entre niveis nao consecutivos no arquivo de distancias: "
				 << lu << " -> " << lv << endl;
			return false;
		}

		if (iu >= FWD[lu].size() || iv >= BWD[lv].size()) {
			cerr << " Id fora da faixa no arquivo de distancias: (" << lu << ","
				 << iu << ") -> (" << lv << "," << iv << ")."
				 << " O arquivo de distancias corresponde a esta instancia?" << endl;
			return false;
		}

		FWD[lu][iu].push_back(make_pair(iv, d));
		BWD[lv][iv].push_back(make_pair(iu, d));

		read++;

	}

	if (declared != read) {
		cerr << " AVISO: o arquivo de distancias declara " << declared
			 << " arcos, mas contem " << read << endl;
	}

	n_arcs = read;

	return true;

}
