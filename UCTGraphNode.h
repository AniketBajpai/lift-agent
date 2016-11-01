#ifndef LIFT_AGENT_UCTGRAPHNODE_H
#define LIFT_AGENT_UCTGRAPHNODE_H

#include <unordered_map>
#include <cmath>
#include "State.h"


class UCTGraphNode {
public:
	State *state;
	int numSimulations;
	unordered_map<int, int> simulations;
	unordered_map<int, double> actionCosts;

	UCTGraphNode(State *state) {
		this->state = new State(*state);    // Create copy of passed state
		numSimulations = 0;
		simulations = unordered_map<int, int>();
		actionCosts = unordered_map<int, double>();
	}

	static int getIntFromAction(Action action) {
		int ae1 = action.first;
		int ae2 = action.second;
		const int M = 10;
		return (ae1 * M + ae2);
	}

	// get action according to UCB policy
	Action getUCBAction() {
		vector<Action> actions = state->getSimpleActions();
		double maxValue = INT_MIN;
		Action bestAction;
		double value;
		for (auto action: actions) {
			const int actionHash = getIntFromAction(action);
			value = actionCosts[actionHash] + C_policy * sqrt(log(simulations[actionHash]) / numSimulations);
			if (maxValue < value) {
				maxValue = value;
				bestAction = action;
			}
		}
		return bestAction;
	}

	// Update node cost based on cost for given action
	void updateCost(Action action, double cost) {
		numSimulations++;
		const int actionHash = getIntFromAction(action);
		double actionCost = actionCosts[actionHash];
		int actionSimulations = simulations[actionHash];
		actionSimulations++;
		double costIncrement = (cost - actionCost) / (double) actionSimulations;
		actionCost += costIncrement;
		simulations[actionHash] = actionSimulations;
		actionCosts[actionHash] = actionCost;
	}
};

#endif //LIFT_AGENT_UCTGRAPHNODE_H
