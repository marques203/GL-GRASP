/*
 * EliteSet.h
 *
 *  Created on: 12 giu 2017
 *      Author: antonio
 */

#include "HDAG.h"

#ifndef ELITESET_H_
#define ELITESET_H_

using namespace std;

class Elite_Set {
public:
	Elite_Set(unsigned D, double th);
	virtual ~Elite_Set();

	unsigned getBestElite() const {
		return best_Elite;
	}

	void setBestElite(unsigned bestElite) {
		best_Elite = bestElite;
	}

	const vector<HDAG*>& getElite() const {
		return Elite;
	}

	const vector<double>& getEliteCost() const {
		return Elite_cost;
	}

	void setElite(const vector<HDAG*>& elite) {
		Elite = elite;
	}

	unsigned getSizeElite() const {
		return size_Elite;
	}

	void setSizeElite(unsigned sizeElite) {
		size_Elite = sizeElite;
	}

	unsigned getWorstElite() const {
		return worst_Elite;
	}

	void setWorstElite(unsigned worstElite) {
		worst_Elite = worstElite;
	}

	bool is_elegible(HDAG &S, double cost);

	double get_diversity_degree(HDAG &S1, HDAG &S2);

	bool add_to_Elite(HDAG &S, double cost);

private:

	double threshold;
	unsigned size_Elite;
	unsigned best_Elite;
	unsigned worst_Elite;
	vector<HDAG *> Elite;
	vector<double> Elite_cost;

};

#endif /* ELITESET_H_ */
