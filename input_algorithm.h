/*
 * input_algorithm.h
 *
 *  Created on: 02 apr 2017
 *      Author: antonio
 */

#include <iostream>
#include <sstream>
#include <iomanip>
#include <type_traits>
#include <string.h>
#include <map>

#include <cstdint>

#ifndef INPUT_ALGORITHM_H_
#define INPUT_ALGORITHM_H_

using namespace std;

unsigned check_algorithm_name(int argc, char *argv[]) {

	char algorithm[20];
	strcpy(algorithm, argv[2]);

	if (strcmp(algorithm, "cplex") == 0 ) return 0;

	if (strcmp(algorithm, "grasp1") == 0 ) return 1;

	if (strcmp(algorithm, "grasp2") == 0 ) return 2;

	if (strcmp(algorithm, "grasp3") == 0 ) return 3;

	if (strcmp(algorithm, "tabu") == 0 ) return 4;

	if (strcmp(algorithm, "lsolver") == 0 ) return 5;

	return 6;
}

bool correct_cplex_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " k-value (arg[3]) must be an integer positive value" << endl;
		exit(0);
	}
	// - check the parameter [6], the time-limit
	int tl = atoi(argv[6]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[6]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}

	return true;

}

bool correct_localsolver_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " k-value (arg[3]) must be an integer positive value" << endl;
		exit(0);
	}
	// - check the parameter [6], the time-limit
	int tl = atoi(argv[6]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[6]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}

	return true;

}

bool correct_grasp_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " k-value (arg[3]) must be an integer positive value" << endl;
		exit(0);
	}
	// - check the parameter [4], the alpha value
	double alpha = atof(argv[4]);
	if (alpha < 0 || alpha > 1) {
		cerr << " apha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
		exit(0);
	}
	// - check the parameter [5], local_search
	char *ls = argv[5];
	if (strcmp(ls, "first") != 0 ) {
		if (strcmp(ls, "best") != 0 ) {
			if (strcmp(ls, "no") != 0 ) {
				if (strcmp(ls, "first-prepro") != 0){
					cerr << " Selected local-search (arg[5]) incorrect!" << endl;
					cerr << " - Possible choices {first, first-prepro, best, no}" << endl;
					exit(0);
				}
			}
		}
	}
	// - check the parameter [6], max_iteration
	int maxit = atoi(argv[6]);
	if (maxit < 1 || maxit > 1000000) {
		cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 1000000]" << endl;
		exit(0);
	}
	// - check the parameter [10], the time-limit
	int tl = atoi(argv[10]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}
	// - check the parameters [11], [12] and [13]:
	// - [11]: type of path relinking
	// - [12]: dimension of the elite set
	// - [13]: threshold for DIVERSITY DEGREE
	if (argc > 11) {
		if (argc < 14) {
			cerr << " the parameters [11], [12], and [13] must be :" << endl;
			cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
					"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
			cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
			cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
			exit(0);
		}
		else {
			char *t_PR = argv[11];
			int dim_ES = atoi(argv[12]);
			double th = atof(argv[13]);
			if (strcmp(t_PR, "f") != 0 && strcmp(t_PR, "F") != 0 &&
				strcmp(t_PR, "b") != 0 && strcmp(t_PR, "B") != 0 &&
				strcmp(t_PR, "m") != 0 && strcmp(t_PR, "M") != 0 ) {
				cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
						"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
				exit(0);
			}
			if (dim_ES < 5 || dim_ES > 20) {
				cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
				exit(0);
			}
			if (th < 0 || th > 1) {
				cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
				exit(0);
			}
		}
	}

	return true;

}

bool correct_tabu_parameter(int argc, char *argv[]) {

	// - check the parameter [3], the k value
	int k = atoi(argv[3]);
	if (k < 0) {
		cerr << " k-value (arg[3]) must be an integer positive value" << endl;
		exit(0);
	}
	// - check the parameter [4], the alpha value
	double alpha = atof(argv[4]);
	if (alpha < 0 || alpha > 1) {
		cerr << " apha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
		exit(0);
	}
	// - check the parameter [5], local_search
	char *ls = argv[5];
	if (strcmp(ls, "first") != 0 ) {
		if (strcmp(ls, "best") != 0 ) {
			if (strcmp(ls, "no") != 0 ) {
				if (strcmp(ls, "first-prepro") != 0){
					cerr << " Selected local-search (arg[5]) incorrect!" << endl;
					cerr << " - Possible choices {first, first-prepro, best, no}" << endl;
					exit(0);
				}
			}
		}
	}
	// - check the parameter [6], max_iteration
/*	int maxit = atoi(argv[6]);				// NEW!!!!!!!!!!!!
	if (maxit < 1 || maxit > 100000) {
		cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 1000]" << endl;
		exit(0);
	}*/
	// - check the parameter [10], the time-limit
	int tl = atoi(argv[10]);
	if (tl < 10 || tl > 3600) {
		cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
		exit(0);
	}
	// - check the parameters [11], [12] and [13]:
	// - [11]: type of path relinking
	// - [12]: dimension of the elite set
	// - [13]: threshold for DIVERSITY DEGREE
	if (argc > 11) {
		if (argc < 14) {
			cerr << " the parameters [11], [12], and [13] must be :" << endl;
			cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
					"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
			cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
			cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
			exit(0);
		}
		else {
			char *t_PR = argv[11];
			int dim_ES = atoi(argv[12]);
			double th = atof(argv[13]);
			if (strcmp(t_PR, "f") != 0 && strcmp(t_PR, "F") != 0 &&
				strcmp(t_PR, "b") != 0 && strcmp(t_PR, "B") != 0 &&
				strcmp(t_PR, "m") != 0 && strcmp(t_PR, "M") != 0 ) {
				cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
						"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
				exit(0);
			}
			if (dim_ES < 5 || dim_ES > 20) {
				cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
				exit(0);
			}
			if (th < 0 || th > 1) {
				cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
				exit(0);
			}
		}
	}

	return true;

}

unsigned input_manager(char argc, char *argv[]) {

	unsigned algorithm = check_algorithm_name(argc, argv);

	switch(algorithm) {
	case(0):
	// "Execute cplex";
	{

		if (correct_cplex_parameter(argc, argv))
			return 0;

 	}
	break;
	case(1):
	{

		if (correct_grasp_parameter(argc, argv))
			return 1;

	}
	break;
	case(2):
	{

		if (correct_grasp_parameter(argc, argv))
			return 2;

	}
	break;
	case(3):
	{

		if (correct_grasp_parameter(argc, argv))
			return 3;

	}
	break;
	case(4):
	{
		if (correct_tabu_parameter(argc, argv))
			return 4;
	}
	break;
	case(5):
	{
		if (correct_localsolver_parameter(argc, argv))
			return 5;
	}
	break;
	}

	return 0;


/*	// Algorithm selected
	char algorithm[20];
	strcpy(algorithm, argv[2]);

	if (strcmp(algorithm, "cplex") != 0 ) {
		if (strcmp(algorithm, "grasp1") != 0 ) {
			if (strcmp(algorithm, "grasp2") != 0 ) {
				cerr << " Selected algorithm (arg[2]) incorrect!" << endl;
				cerr << " - Possible choices {cplex, grasp1, grasp2}" << endl;
				exit(0);
			}
			else {
				// Name of the algorithm is "grasp2"
				// - check the parameter [3], the k value
				int k = atoi(argv[3]);
				if (k < 0) {
					cerr << " k-value (arg[3]) must be an integer positive value" << endl;
					exit(0);
				}
				// - check the parameter [4], the alpha value
				double alpha = atof(argv[4]);
				if (alpha < 0 || alpha > 1) {
					cerr << " apha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
					exit(0);
				}
				// - check the parameter [5], local_search
				char *ls = argv[5];
				if (strcmp(ls, "first") != 0 ) {
					if (strcmp(ls, "best") != 0 ) {
						if (strcmp(ls, "no") != 0 ) {
							if (strcmp(ls, "first-prepro") != 0){
								cerr << " Selected local-search (arg[5]) incorrect!" << endl;
								cerr << " - Possible choices {first, first-prepro, best, no}" << endl;
								exit(0);
							}
						}
					}
				}
				// - check the parameter [6], max_iteration
				int maxit = atoi(argv[6]);
				if (maxit < 1 || maxit > 100000) {
					cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 1000]" << endl;
					exit(0);
				}
				// - check the parameter [10], the time-limit
				int tl = atoi(argv[10]);
				if (tl < 10 || tl > 3600) {
					cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
					exit(0);
				}
				// - check the parameters [11], [12] and [13]:
				// - [11]: type of path relinking
				// - [12]: dimension of the elite set
				// - [13]: threshold for DIVERSITY DEGREE
				if (argc > 11) {
					if (argc < 14) {
						cerr << " the parameters [11], [12], and [13] must be :" << endl;
						cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
								"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
						cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
						cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
						exit(0);
					}
					else {
						char *t_PR = argv[11];
						int dim_ES = atoi(argv[12]);
						double th = atof(argv[13]);
						if (strcmp(t_PR, "f") != 0 && strcmp(t_PR, "F") != 0 &&
							strcmp(t_PR, "b") != 0 && strcmp(t_PR, "B") != 0 &&
							strcmp(t_PR, "m") != 0 && strcmp(t_PR, "M") != 0 ) {
							cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
									"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
							exit(0);
						}
						if (dim_ES < 5 || dim_ES > 20) {
							cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
							exit(0);
						}
						if (th < 0 || th > 1) {
							cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
							exit(0);
						}
					}
				}
				return 2;
			}
		}
		else {
			// Name of the algorithm is "grasp1"
			// - check the parameter [3], the k value
			int k = atoi(argv[3]);
			if (k < 0) {
				cerr << " k-value (arg[3]) must be an integer positive value" << endl;
				exit(0);
			}
			// - check the parameter [4], the alpha value
			double alpha = atof(argv[4]);
			if (alpha < 0 || alpha > 1) {
				cerr << " apha-value (arg[4]) must be a real nonnegative value in the range [0,1]" << endl;
				exit(0);
			}
			// - check the parameter [5], local_search
			char *ls = argv[5];
			if (strcmp(ls, "first") != 0 ) {
				if (strcmp(ls, "best") != 0 ) {
					if (strcmp(ls, "no") != 0 ) {
						if (strcmp(ls, "first-prepro") != 0){
							cerr << " Selected local-search (arg[5]) incorrect!" << endl;
							cerr << " - Possible choices {first, first-prepro, best, no}" << endl;
							exit(0);
						}
					}
				}
			}
			// - check the parameter [6], max_iteration
			int maxit = atoi(argv[6]);
			if (maxit < 1 || maxit > 100000) {
				cerr << " max-it (arg[6]) must be an integer positive value in the range [10, 1000]" << endl;
				exit(0);
			}
			// - check the parameter [10], the time-limit
			int tl = atoi(argv[10]);
			if (tl < 10 || tl > 3600) {
				cerr << " time-limit (arg[10]) must be an integer positive value in the range [10, 3600]" << endl;
				exit(0);
			}
			// - check the parameters [11], [12] and [13]:
			// - [11]: type of path relinking
			// - [12]: dimension of the elite set
			// - [13]: threshold for DIVERSITY DEGREE
			if (argc > 11) {
				if (argc < 14) {
					cerr << " the parameters [11], [12], and [13] must be :" << endl;
					cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
							"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
					cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
					cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
					exit(0);
				}
				else {
					char *t_PR = argv[11];
					int dim_ES = atoi(argv[12]);
					double th = atof(argv[13]);
					if (strcmp(t_PR, "f") != 0 && strcmp(t_PR, "F") != 0 &&
						strcmp(t_PR, "b") != 0 && strcmp(t_PR, "B") != 0 &&
						strcmp(t_PR, "m") != 0 && strcmp(t_PR, "M") != 0 ) {
						cerr << " - [11]: type of path relinking must be a char in {f(or F), b(or B), m(or M)} "
								"(f/F: forward, b/B: backward, m/M: mixed)" << endl;
						exit(0);
					}
					if (dim_ES < 5 || dim_ES > 20) {
						cerr << " - [12]: dimension of the elite set, an integer positive value in the range [5, 20]" << endl;
						exit(0);
					}
					if (th < 0 || th > 1) {
						cerr << " - [13]: threshold for access to the elite set, a real value in the range [0, 1]" << endl;
						exit(0);
					}
				}
			}
			return 1;
		}
	}
	else {
		// Name of the algorithm is "cplex"
		// - check the parameter [3], the k value
		int k = atoi(argv[3]);
		if (k < 0) {
			cerr << " k-value (arg[3]) must be an integer positive value" << endl;
			exit(0);
		}
		// - check the parameter [6], the time-limit
		int tl = atoi(argv[6]);
		if (tl < 10 || tl > 3600) {
			cerr << " time-limit (arg[6]) must be an integer positive value in the range [10, 3600]" << endl;
			exit(0);
		}
		return 0;
	}

	return 0;
	*/

}


#endif /* INPUT_ALGORITHM_H_ */
