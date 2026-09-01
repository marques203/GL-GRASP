/*
 * EliteSet.cpp
 *
 *  Created on: 12 giu 2017
 *      Author: antonio
 */

#include "EliteSet.h"

Elite_Set::Elite_Set(unsigned D, double th) : threshold(0), size_Elite(10), best_Elite(0), worst_Elite(0) {
	// TODO Auto-generated constructor stub

	size_Elite = D;
	threshold = th;

}

bool Elite_Set::add_to_Elite(HDAG &S, double cost) {

	if (Elite.size() == 0) {

		Elite.push_back(&S);
		Elite_cost.push_back(cost);
		best_Elite = 0;

		return true;

	}

	if (S.get_cost() < Elite_cost[best_Elite] || is_elegible(S, cost)) {
		if (Elite.size() == size_Elite) {
			// Replace the worst element in Elite
			HDAG &tmp = *Elite[worst_Elite];
			Elite[worst_Elite] = NULL;
			Elite[worst_Elite] = &S;
			Elite_cost[worst_Elite] = cost;
			if (Elite_cost[worst_Elite] < Elite_cost[best_Elite]) {
				best_Elite = worst_Elite;
			}
			double worst = 0;
			// find new worst
			for (unsigned i = 0; i < Elite.size(); ++i) {
				if (Elite_cost[i] > worst) {
					worst = Elite_cost[i];
					worst_Elite = i;
				}
			}
			free(&tmp);
		}
		else {
			Elite.push_back(&S);
			Elite_cost.push_back(cost);
			if (cost < Elite_cost[best_Elite]) {
				best_Elite = Elite.size() - 1;
			}
		}
		return true;
	}

	return false;

}

bool Elite_Set::is_elegible(HDAG &S, double cost) {

	for (unsigned i = 0; i < Elite.size(); ++i) {
		if (get_diversity_degree(*Elite[i], S) >= threshold) {
			if (Elite.size() == size_Elite && cost >= Elite_cost[worst_Elite]) {
				return false;
			}
		}
		else {
			return false;
		}
	}

	return true;

}


double Elite_Set::get_diversity_degree(HDAG &S1, HDAG &S2) {

	auto &LEVELS = S1.getLEVELS();
	auto &IDs_S1 = S1.getIDs();
	auto &IDs_S2 = S2.getIDs();
	unsigned diversity = 0;

	for (unsigned l = 0; l < LEVELS.size(); ++l) {
		for (unsigned i = 0; i < LEVELS[l].size(); ++i) {
			if (IDs_S1[l][i] != IDs_S2[l][i])
				diversity++;
		}
	}

	return (double)((double)diversity / (double)S1.get_total_nodes());

}


Elite_Set::~Elite_Set() {
	// TODO Auto-generated destructor stub
}

