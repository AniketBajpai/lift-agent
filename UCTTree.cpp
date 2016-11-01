#include "UCTTree.h"

UCTTree::UCTTree(State *startState) {
	this->startState = startState;
	this->rootNode = new UCTNode();
	this->rootNode->parent = nullptr;
}

// select node according to tree policy
vector<Action>& UCTTree::getActionSequence() {
	vector<Action> actions;
	UCTNode *chosenNode = this->rootNode;
	while (chosenNode->children.size() > 0) {
		auto next = chosenNode->chooseNext();
		Action action = next.first;
		action.push_back(action);
		chosenNode = next.second;
	}
	return actions;
}

UCTNode* UCTTree::getTreeNode(vector<Action>& actions) {
	UCTNode* currNode = this->rootNode;
	for(auto action: actions) {
		currNode = currNode->getChildNode(action);
	}
	return currNode;
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
	UCTNode *newNode = new UCTNode();
	newNode->parent = node;
	node->children.push_back(make_pair(newAction, newNode));
	return newNode;
}

// run monte carlo simulation from chosen node and backprop the value in tree
void UCTTree::simulate(UCTNode *node) {
	State *currentState = node->state;
	double simValue;
	for (int i = 0; i < NUM_SIM; ++i) {
		double simValueCurr = runSimulation(currentState,  ,SIM_STATE_EPOCHS);
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
		vector<Action> actions = this->select();
		UCTNode *newNode = this->expand(selectedNode);
		this->simulate(newNode);
	}
}

double UCTTree::runSimulation(State* startState, vector<Action>& actions, int epochs) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));

	for(auto action: actions) {
		startState->simulateStep(num_out);
		cost += startState->applyAction(action, num_out);   // Initial action cost
		cost += startState->getWaitCost(num_out);
	}

	ElevatorState  prevState, currState;
	while (epochs > 0) {
		prevState = startState->elevator1.elevatorState;
		Action chosenAction = startState->getPolicyAction();
		startState->simulateStep(num_out);  // simulator step
		currState = startState->elevator1.elevatorState;
		cost += startState->applyAction(chosenAction, num_out);  // change lift state according to action, add action cost
		cost += startState->getWaitCost(num_out);  // Add wait cost
		bool isNewEpoch = (prevState == FULL and currState == EMPTY);
		if(isNewEpoch)
			epochs--;
	}

	return cost;
}

