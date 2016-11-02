#include "UCTGraph.h"

vector<UCTGraphNode *> UCTGraph::nodeList = vector<UCTGraphNode *>();

UCTGraphNode *UCTGraph::getNodeFromState(State *state) {
	for (auto node: nodeList) {
		if (node->state == state) {
			return node;
		}
	}
	return nullptr;
}

void UCTGraph::run(State *startState, int depth) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));

	runhelper(startState, depth, num_out);
}

double UCTGraph::runhelper(State *state, int depth, int num_out[N + 1][N + 1]) {
	if (depth == 0) {
		return 0;
	}
	UCTGraphNode *node = getNodeFromState(state);
	if (node) {      // State found
		double cost = 0;
		// Get action according to (UCB) policy
		Action ucbAction = node->getUCBAction();

		// Simulate action to change state
		// TODO: check order of adding costs
		cost += state->applyAction(ucbAction, num_out); // add action cost
		cost += state->getWaitCost(num_out);  // Add wait cost

		// Call runhelper on new state to get cost
		cost += runhelper(state, depth - 1, num_out);

		// Update cost and number of simulations
		node->updateCost(ucbAction, cost);

		return cost;
	}
	else {      // State not found in list
		UCTGraphNode *newNode = new UCTGraphNode(state);

		// get vector of actions from state;
		vector<Action> actions = state->getSimpleActions();

		// initialize number of simulations and costs for actions
		for (auto action:actions) {
			const int actionHash = UCTGraphNode::getIntFromAction(action);
			newNode->simulations[actionHash] = 0;
			newNode->actionCosts[actionHash] = 0;
		}

		// run simulation for action from state
		for (auto action:actions) {
			const int actionHash = UCTGraphNode::getIntFromAction(action);
			newNode->simulations[actionHash] = 1;
			newNode->actionCosts[actionHash] = runSimulation(state, action, depth);   // Update action costs
		}
		nodeList.push_back(newNode);
	}
}

double UCTGraph::runSimulation(State *startState, Action action, int depth) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));

	State *state = new State(*startState);
	state->simulateStep(num_out);
	cost += state->applyAction(action, num_out);   // Initial action cost
	cost += state->getWaitCost(num_out);

	while (depth--) {
		Action chosenAction = state->getPolicyAction().first;
		state->simulateStep(num_out);  // simulator step
		cost += state->applyAction(chosenAction, num_out);  // change lift state according to action, add action cost
		cost += state->getWaitCost(num_out);  // Add wait cost
	}
	cost += state->getPolicyAction().second;

	return cost;
}

vector<pair<Action, double> > UCTGraph::getBaseCosts(State *state, int depth) {
	vector<pair<Action, double> > baseCosts;
	vector<Action> actions = state->getSimpleActions();
	double actionCost;
	for (auto action: actions) {
		actionCost = runSimulation(state, action, depth);
		baseCosts.push_back(make_pair(action, actionCost));
	}

	return baseCosts;
}




