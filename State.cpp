#include <vector>
#include <cstdlib>
#include <cstring>
#include <random>
#include "State.h"
#include "utility.h"

vector<Action> State::getActions() {
	vector<ElevatorAction> actions1 = elevator1.getActions();
	vector<ElevatorAction> actions2 = elevator2.getActions();
	vector<Action> actions;
	for (const auto &&action1: actions1) {
		for (const auto &&action2: actions2) {
			actions.push_back(make_pair(action1, action2));
		}
	}
	return actions;
}

State *State::getResState(Action action) {
	State *resState = new State(*this);

	return resState;
}


pair<Action, State *> State::getRandomNextState() {
	vector<Action> actions = getActions();
	int l = actions.size();
	int i = rand() % l;
	return make_pair(actions[i], getResState(actions[i]));
}

// get action according to policy
Action State::getPolicyAction() {

	return std::pair<ElevatorAction, ElevatorAction>();
}

void State::simulateStep(int num_out[N + 1][N + 1]) {
	// simulate entry of people in system
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);
	double x = distribution(generator);

	int i;  // stores floor number
	int j;  // stores alight floor number

	if (x < p) {   // Person arrives with probability p
		// Decide which floor person arrives
		double y = distribution(generator);
		double delta = (1.0 - q) / (double) (N - 1);
		if (y < q) {
			i = 1;
		}
		else {
			i = (int) ((y - q) / delta) + 2;
		}

		// Decide floor number where person alights
		double z = distribution(generator);
		if (i == 1) {
			delta = 1.0 / (double) (N - 1);
			j = (int) (z / delta) + 2;
		}
		else {
			if (z < r) {
				j = 1;
			}
			else {
				delta = (1.0 - r) / (double) (N - 2);
				j = (int) (z / delta) + 2;
				if (j >= i)
					j++;
			}
		}
	}

	num_out[i][j]++;
}

double State::applyAction(Action action, int num_out[N + 1][N + 1]) {
	double cost = 0;
	cost += this->elevator1.applyAction(action.first, num_out);
	cost += this->elevator1.applyAction(action.second, num_out);
	return cost;
}


double State::getWaitCost(Action action, int num_out[N + 1][N + 1]) {
	double cost = 0;
	cost += (this->elevator1.getNumberOfPeople() * WAIT_COST);
	cost += (this->elevator2.getNumberOfPeople() * WAIT_COST);
	return cost;
}


double State::runSimulation(int epochs) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));

	State *startState = new State(*this);

	Action action = startState->getPolicyAction();
	startState->simulateStep(num_out);  // simulator step
	startState->applyAction(action, num_out);    // change lift state according to action, add action cost

	cost += startState->getWaitCost(action, num_out);  // Add wait cost


	return cost;
}

