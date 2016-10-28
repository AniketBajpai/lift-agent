#ifndef LIFT_AGENT_UCTNODE_H
#define LIFT_AGENT_UCTNODE_H

#include "State.h"

class UCTNode {
public:
	UCTNode *parent;
	vector<pair<Action, UCTNode *> > children;
	State *state;
	int num_simulations;
	double value;

	UCTNode(State *state);

	double getTreeValue();

	void updateValue(double value);

	pair<Action, UCTNode *> chooseNext();
};

#endif //LIFT_AGENT_UCTNODE_H
