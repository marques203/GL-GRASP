/*
 * test_glgrasp_units.cpp
 *
 * Testes de unidade das duas pecas de logica introduzidas no GL-GRASP para
 * alinhar a implementacao ao codigo de referencia (Cb_*.py) de Charytitsch e
 * Nascimento (2026):
 *
 *   - gl_rcl_pick_uniform : sorteio UNIFORME sobre a RCL achatada (item 1.1)
 *   - gl_func_ale         : alvo uniforme entre {min, max, ponto medio} do
 *                           Step 5 / Figura 3 (item 1.2)
 *
 * As duas sao funcoes livres declaradas em GLGRASP.h, logo os testes abaixo
 * exercitam exatamente o codigo que roda em producao.
 *
 * Compilar e rodar:
 *     gmake test
 */

#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <cmath>

// Inclui o header DE PRODUCAO: gl_func_ale (inline) e a declaracao de
// gl_rcl_pick_uniform vem dali, nao de uma copia. Se a assinatura mudar, este
// arquivo deixa de compilar ou de linkar, o que tambem e informacao.
#include "../GLGRASP.h"

using namespace std;

static int failures = 0;
static int checks = 0;

static void check(bool ok, const char *name, const char *detail = "") {
	checks++;
	if (ok) {
		printf("  [ ok ] %s %s\n", name, detail);
	} else {
		printf("  [FAIL] %s %s\n", name, detail);
		failures++;
	}
}

// RCL de teste com niveis de tamanhos bem desiguais: e exatamente nessa
// situacao que o sorteio "nivel primeiro" da classe base se afasta do
// uniforme. Total = 1 + 1 + 10 + 0 + 4 = 16 entradas.
static const unsigned LEVEL_SIZES[] = { 1, 1, 10, 0, 4 };
static const unsigned N_LEVELS = 5;
static const unsigned TOTAL = 16;

static vector<vector<pair<unsigned, unsigned> > > make_rcl() {
	vector<vector<pair<unsigned, unsigned> > > RCL(N_LEVELS);
	for (unsigned l = 0; l < N_LEVELS; ++l)
		for (unsigned j = 0; j < LEVEL_SIZES[l]; ++j)
			RCL[l].push_back(make_pair(l, j));
	return RCL;
}

// ---------------------------------------------------------------------------
// 1.1 (a) Dreno: esvaziar a RCL tem de devolver cada entrada EXATAMENTE uma
// vez, sem indice fora de faixa e sem repeticao.
// ---------------------------------------------------------------------------
static void test_pick_drains_exactly() {

	bool all_ok = true;

	for (unsigned seed = 1; seed <= 200 && all_ok; ++seed) {

		MTRand rng(seed);
		vector<vector<pair<unsigned, unsigned> > > RCL = make_rcl();

		multiset<pair<unsigned, unsigned> > seen;

		for (unsigned t = TOTAL; t > 0; --t) {

			pair<unsigned, unsigned> p = gl_rcl_pick_uniform(RCL, t, rng);

			if (p.first >= N_LEVELS || p.second >= LEVEL_SIZES[p.first]) {
				all_ok = false;
				printf("        seed %u: par fora de faixa (%u,%u)\n",
						seed, p.first, p.second);
				break;
			}

			seen.insert(p);

		}

		unsigned remaining = 0;
		for (unsigned l = 0; l < N_LEVELS; ++l)
			remaining += RCL[l].size();

		if (remaining != 0) {
			all_ok = false;
			printf("        seed %u: sobraram %u entradas\n", seed, remaining);
		}

		if (seen.size() != TOTAL) {
			all_ok = false;
			printf("        seed %u: %u extracoes em vez de %u\n",
					seed, (unsigned) seen.size(), TOTAL);
		}

		// cada entrada uma unica vez
		for (unsigned l = 0; l < N_LEVELS && all_ok; ++l)
			for (unsigned j = 0; j < LEVEL_SIZES[l]; ++j)
				if (seen.count(make_pair(l, j)) != 1) {
					all_ok = false;
					printf("        seed %u: (%u,%u) saiu %u vezes\n", seed, l, j,
							(unsigned) seen.count(make_pair(l, j)));
				}

	}

	check(all_ok, "1.1(a) dreno da RCL devolve cada entrada exatamente uma vez",
			"(200 sementes)");

}

// ---------------------------------------------------------------------------
// 1.1 (b) Uniformidade: a frequencia de cada uma das 16 entradas tem de ficar
// proxima de 1/16. Para contraste, o mesmo teste e aplicado ao esquema
// "nivel primeiro" da classe base, que deve FALHAR.
// ---------------------------------------------------------------------------
static pair<unsigned, unsigned> pick_level_first(
		vector<vector<pair<unsigned, unsigned> > > &RCL, MTRand &rng) {

	// Replica de GRASP::extract_from_RCL, so para medir o contraste.
	unsigned l, i;
	while (1) {
		l = rng.randInt(RCL.size() - 1);
		if (RCL[l].size() > 0)
			break;
	}
	i = rng.randInt(RCL[l].size() - 1);
	return RCL[l][i];

}

static void test_pick_is_uniform() {

	const unsigned N = 1600000;   // 100000 esperados por entrada

	map<pair<unsigned, unsigned>, unsigned> freq_uniform;
	map<pair<unsigned, unsigned>, unsigned> freq_level_first;

	MTRand rng_a(12345);
	MTRand rng_b(12345);

	for (unsigned t = 0; t < N; ++t) {

		vector<vector<pair<unsigned, unsigned> > > RCL = make_rcl();
		freq_uniform[gl_rcl_pick_uniform(RCL, TOTAL, rng_a)]++;

		vector<vector<pair<unsigned, unsigned> > > RCL2 = make_rcl();
		freq_level_first[pick_level_first(RCL2, rng_b)]++;

	}

	double expected = (double) N / (double) TOTAL;

	double worst_uniform = 0.0, worst_level_first = 0.0;

	for (unsigned l = 0; l < N_LEVELS; ++l) {
		for (unsigned j = 0; j < LEVEL_SIZES[l]; ++j) {

			pair<unsigned, unsigned> key = make_pair(l, j);

			double du = fabs(freq_uniform[key] - expected) / expected;
			if (du > worst_uniform) worst_uniform = du;

			double dl = fabs(freq_level_first[key] - expected) / expected;
			if (dl > worst_level_first) worst_level_first = dl;

		}
	}

	char buf[256];

	// 1.6e6 amostras, p = 1/16: sigma/esperado ~ 0.3%. 2% e ~6 sigma.
	snprintf(buf, sizeof(buf), "(desvio maximo %.2f%%, limite 2%%)",
			100.0 * worst_uniform);
	check(worst_uniform < 0.02, "1.1(b) sorteio uniforme sobre a RCL inteira", buf);

	// Confirma que o teste tem poder: o esquema antigo erra por muito.
	snprintf(buf, sizeof(buf), "(desvio maximo %.0f%% no esquema nivel-primeiro)",
			100.0 * worst_level_first);
	check(worst_level_first > 0.50,
			"1.1(b) contraste: esquema da classe base NAO e uniforme", buf);

	// Nenhuma entrada pode ficar sem ser sorteada.
	bool all_hit = true;
	for (unsigned l = 0; l < N_LEVELS; ++l)
		for (unsigned j = 0; j < LEVEL_SIZES[l]; ++j)
			if (freq_uniform[make_pair(l, j)] == 0)
				all_hit = false;
	check(all_hit, "1.1(b) toda entrada da RCL tem probabilidade positiva", "");

}

// ---------------------------------------------------------------------------
// 1.1 (c) Determinismo: mesma semente, mesma sequencia.
// ---------------------------------------------------------------------------
static void test_pick_deterministic() {

	vector<pair<unsigned, unsigned> > a, b;

	for (int run = 0; run < 2; ++run) {

		MTRand rng(777);
		vector<vector<pair<unsigned, unsigned> > > RCL = make_rcl();

		for (unsigned t = TOTAL; t > 0; --t)
			(run == 0 ? a : b).push_back(gl_rcl_pick_uniform(RCL, t, rng));

	}

	check(a == b, "1.1(c) sorteio reprodutivel com semente fixa", "");

}

// ---------------------------------------------------------------------------
// 1.2 (a) Suporte do func_ale: exatamente {a, b, (a+b)/2}.
// ---------------------------------------------------------------------------
static void test_func_ale_support() {

	MTRand rng(42);

	unsigned a = 4, b = 11;
	unsigned mid = (a + b) / 2;   // 7

	set<unsigned> seen;
	for (unsigned t = 0; t < 20000; ++t)
		seen.insert(gl_func_ale(a, b, rng));

	set<unsigned> expected;
	expected.insert(a);
	expected.insert(b);
	expected.insert(mid);

	check(seen == expected,
			"1.2(a) func_ale gera exatamente {min, max, ponto medio}", "");

	// O ponto medio DENTRO de um nivel e justamente o que a leitura do texto
	// (POS_RULE_TEXT) nao produz; aqui ele tem de aparecer.
	check(seen.count(mid) == 1, "1.2(a) ponto medio pertence ao suporte", "");

}

// ---------------------------------------------------------------------------
// 1.2 (b) Uniformidade do func_ale entre as tres opcoes.
// ---------------------------------------------------------------------------
static void test_func_ale_uniform() {

	MTRand rng(2024);

	const unsigned N = 900000;
	unsigned a = 2, b = 9;              // a, b e (a+b)/2 = 5 sao distintos
	map<unsigned, unsigned> freq;

	for (unsigned t = 0; t < N; ++t)
		freq[gl_func_ale(a, b, rng)]++;

	double expected = (double) N / 3.0;
	double worst = 0.0;
	unsigned opts[3] = { a, b, (a + b) / 2 };

	for (int j = 0; j < 3; ++j) {
		double d = fabs(freq[opts[j]] - expected) / expected;
		if (d > worst) worst = d;
	}

	char buf[128];
	snprintf(buf, sizeof(buf), "(desvio maximo %.2f%%, limite 1%%)", 100.0 * worst);
	check(worst < 0.01, "1.2(b) func_ale uniforme entre as tres opcoes", buf);

}

// ---------------------------------------------------------------------------
// 1.2 (c) Casos de borda do func_ale.
// ---------------------------------------------------------------------------
static void test_func_ale_edges() {

	MTRand rng(9);

	// a == b: unico valor possivel
	bool ok_eq = true;
	for (unsigned t = 0; t < 1000; ++t)
		if (gl_func_ale(7, 7, rng) != 7)
			ok_eq = false;
	check(ok_eq, "1.2(c) func_ale com a == b devolve sempre a", "");

	// Posicao do vizinho mais proximo pode ser MAIOR que a do mais distante:
	// o ponto medio tem de ser o mesmo nos dois sentidos.
	set<unsigned> s1, s2;
	for (unsigned t = 0; t < 20000; ++t) {
		s1.insert(gl_func_ale(3, 10, rng));
		s2.insert(gl_func_ale(10, 3, rng));
	}
	check(s1 == s2, "1.2(c) func_ale e simetrico em (a, b)", "");

	// Sem overflow em posicoes grandes (a soma e feita em unsigned long).
	unsigned big_a = 4000000000u, big_b = 4000000002u;
	bool ok_big = true;
	for (unsigned t = 0; t < 1000; ++t) {
		unsigned r = gl_func_ale(big_a, big_b, rng);
		if (r != big_a && r != big_b && r != 4000000001u)
			ok_big = false;
	}
	check(ok_big, "1.2(c) func_ale sem overflow em posicoes grandes", "");

}

int main() {

	printf("\n=== Testes de unidade: GL-GRASP, itens 1.1 e 1.2 ===\n\n");

	test_pick_drains_exactly();
	test_pick_is_uniform();
	test_pick_deterministic();
	test_func_ale_support();
	test_func_ale_uniform();
	test_func_ale_edges();

	printf("\n%d verificacoes, %d falha(s)\n\n", checks, failures);

	return failures == 0 ? 0 : 1;

}
