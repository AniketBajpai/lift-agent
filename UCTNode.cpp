#include <cmath>
#include "UCTNode.h"

UCTNode::UCTNode() {
	this->num_simulations = 0;
	this->value = 0;
}

UCTNode* UCTNode::getChildNode(Action action) {
	for(auto child: this->children) {
		if(child.first == action)
			return child.second;
	}
}

void UCTNode::updateValue(double value) {
	if (num_simulations == 0) {
		this->value = value;
	}
	else {
		num_simulations++;
		this->value += (value / num_simulations);
	}
}

double UCTNode::getTreeValue() {
	assert(this->parent != nullptr);
	return (value + C_policy * sqrt(log(num_simulations) / parent->num_simulations));
}

pair<Action, UCTNode *> UCTNode::chooseNext() {
	pair<Action, UCTNode *> chosen;
	double maxValue = INT_MIN;

	assert(children.size() > 0);
	for (auto &&next: children) {
		double nextValue = next.second->getTreeValue();
		if (nextValue > maxValue) {
			maxValue = nextValue;
			chosen = next;
		}
	}
	return chosen;
}

