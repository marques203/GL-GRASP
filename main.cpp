#include <string>
#include <limits>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "input_algorithm.h"
#include "lib/MTRand.h"
#include "HDAG.h"
#ifdef WITH_CPLEX
#include "Cplex.h"
#endif
#include "GRASPv1.h"
#include "GRASPv2.h"
#include "GRASPv3.h"
#include "GLGRASP.h"
#include "TABU.h"
#include <time.h>
#include <ios>
#include "chrono"
#ifdef WITH_LOCALSOLVER
#include "Localsolver.h"
#endif
#include <sstream>


using namespace std;

int main(int argc, char *argv[]) {

	unsigned input = input_manager(argc, argv);
	unsigned k = atoi(argv[3]);

	HDAG I; // Instance of the problem
	I.read_instance(argv[1]); // read the instance

//	I.read_jesusFormat(argv[1]);

//	I.printHDAG(I, "input.tex");
//	I.printHDAGgraphviz(I, "input.txt");

//	cerr << " k = " << argv[3] << " - inst = " << argv[7] << endl;

	switch (input) {
	case(0):
	// "Execute cplex";
	{
#ifdef WITH_CPLEX
		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {
			Cplex C(argv);
			C.cplex(I, argc, argv);
		}
#else
		cerr << "This build was compiled without CPLEX support (WITH_CPLEX not defined)." << endl;
#endif
 	}
	break;
	case(1):
	{

		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {
			GRASP_v1 G(I, argv);
			G.algorithm(I, argc, argv);
		}
	}
	break;
	case(2):
	{

		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {
			GRASP_v2 G(I, argv);
			G.algorithm(I, argc, argv);
		}

	}
	break;
	case(3):
	{

		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {

			GRASP_v3 G(I, argv);
			G.algorithm(I, argc, argv);
//			stringstream outfile;
//			outfile << argv[7] << "_k=" << argv[3] << ".tex";
//			stringstream outfileGraphViz;
//			outfileGraphViz << argv[7] << "_k=" << argv[3] << ".txt";
//			G.get_best_solution().printHDAG(G.get_best_solution(), outfile.str());
//			G.get_best_solution().printHDAGgraphviz(G.get_best_solution(), outfileGraphViz.str());
		}


	}
	break;
	case(4):
	{

		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {
			TABU T(I, argv);
			T.algorithm(I, argc, argv);
		}

	}
	break;
	case(6):
	{

		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {

			GLGRASP G(I, argv);

			// argv[11]: arquivo de distancias gerado pelo modulo Python
			if (G.load_distances(argv[11], I)) {

				// argv[12] (opcional): eta_max, iteracoes sem melhora
				if (argc > 12)
					G.set_eta_max(atoi(argv[12]));

				G.algorithm(I, argc, argv);

			}

		}

	}
	break;
	case(5):
	{
#ifdef WITH_LOCALSOLVER
		if ((int)k <= I.get_MIN_INCREMENTAL_NUMBER(I)) {

			Localsolver LS(argv);
			LS.localsolver(I, argc, argv);
		}
#else
		cerr << "This build was compiled without Localsolver support (WITH_LOCALSOLVER not defined)." << endl;
#endif
	}
	break;
	}


	return 0;

}
