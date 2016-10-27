#include <vector>
#include <cmath>
#include <cassert>
#include "UCTNode.h"

UCTNode::UCTNode(State* state) {
	this->num_simulations = 0;
	vector children();  // empty vector
}

void UCTNode::updateValue(double value) {
	if(num_simulations == 0) {
		this->value = value;
	}
	else {
		this->value += (value/num_simulations);
	}
	num_simulations++;
}

double UCTNode::getTreeValue() {
	return (value + C_policy * sqrt(log(num_simulations)/parent->num_simulations);
}

std::pair<Action, UCTNode *> UCTNode::chooseNext() {
	pair<Action, UCTNode *> chosen;
	double maxValue = INT_MIN;

	assert(children.size() > 0);
	for(auto&& next: children) {
		double nextValue = next.second->getTreeValue();
		if(nextValue > maxValue) {
			maxValue = nextValue;
			chosen = next;
		}
	}
	return chosen;
}

