#ifndef LIFT_AGENT_UCTTREE_H
#define LIFT_AGENT_UCTTREE_H

#include "UCTNode.h"

class UCTTree {
public:
	UCTNode *rootNode;
	State* startState;
	const int NEW_STATE_EPOCHS = 1;
	const int NUM_SIM = 5;
	const int SIM_STATE_EPOCHS = 3;
	const int TREE_ITER = 10;

	UCTTree(State *startState);

	vector<Action>& getActionSequence();

	UCTNode* getTreeNode(vector<Action>& actions);

	UCTNode *select();

	UCTNode *expand(UCTNode *node);

	void simulate(UCTNode *node);

	void run();

	double runSimulation(State* startState, Action action, int epochs);
};

#endif //LIFT_AGENT_UCTTREE_H
