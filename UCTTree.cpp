#include <vector>
#include "UCTTree.h"

UCTTree::UCTTree(State *startState) {
	this->rootNode = new UCTNode(startState);
	this->rootNode->parent = nullptr;
}

// select node according to tree policy
UCTNode *UCTTree::select() {
	UCTNode *chosenNode = this->rootNode;
	while (chosenNode->children.size() > 0) {
		chosenNode = chosenNode->chooseNext().second;
	}
	return chosenNode;
}

// expand chosen node
UCTNode *UCTTree::expand(UCTNode *node) {
	State *currentState = node->state;
	State *newState;
	Action newAction;

//	// greedy expansion
//	double maxValue = INT_MIN;
//	vector<Action> actions = currentState->getActions();
//	for(auto action: actions) {
//		State* resState = currentState->getResState(action);
//		double stateValue = resState->runSimulation(NEW_STATE_EPOCHS);
//		if(stateValue > maxValue) {
//			maxValue = stateValue;
//			newState = resState;
//			newAction = action;
//		}
//	}

	// random expansion
	auto &&random = node->state->getRandomNextState();
	newAction = random.first;
	newState = random.second;

	UCTNode *newNode = new UCTNode(newState);
	newNode->parent = node;
	node->children.push_back(make_pair(newAction, newNode));
	return newNode;
}

// run monte carlo simulation from chosen node and backprop the value in tree
void UCTTree::simulate(UCTNode *node) {
	State *currentState = node->state;
	double simValue;
	for (int i = 0; i < NUM_SIM; ++i) {
		double simValueCurr = currentState->runSimulation(SIM_STATE_EPOCHS);
		simValue += simValueCurr;
	}
	simValue = simValue / NUM_SIM;
	node->updateValue(simValue);

	// backprop
	UCTNode *currNode = node;
	while (currNode->parent) {
		currNode->parent->updateValue(simValue);
		currNode = currNode->parent;
	}
}

void UCTTree::run() {
	for (int i = 0; i < TREE_ITER; ++i) {
		UCTNode *selectedNode = this->select();
		UCTNode *newNode = this->expand(selectedNode);
		this->simulate(newNode);
	}
}

