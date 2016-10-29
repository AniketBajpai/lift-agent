#include <random>
#include "State.h"

queue<ElevatorAction> actionq1;
queue<ElevatorAction> actionq2;


State::State() {
	Elevator elevator1();
	Elevator elevator2();
	memset(time_up, 0, sizeof(time_up));
	memset(time_down, 0, sizeof(time_down));
}


State *State::getResState(Action action) {
	State *resState = new State(*this);

	return resState;
}

//pair<Action, State *> State::getRandomNextState() {
//	vector<Action> actions = getActions();
//	int l = actions.size();
//	int i = rand() % l;
//	return make_pair(actions[i], getResState(actions[i]));
//}

int *State::getDistanceArr(Elevator &elevator, ElevatorAction elevatorAction) {
	int* distance = new int[2 * N + 1];
	int position_new;
	switch (elevatorAction) {
		case AU:
			assert(elevator.position < N);
			position_new = elevator.position + 1;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][1][1];  // Up
				distance[2 * i - 1] = d[position_new][i][1][0];    // Down
			}
			break;
		case AOU:
			assert(elevator.position < N);
			position_new = elevator.position;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][1][1];  // Up
				distance[2 * i - 1] = d[position_new][i][1][0];    // Down
			}
			break;
		case AD:
			assert(elevator.position > 1);
			position_new = elevator.position - 1;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][0][1];  // Up
				distance[2 * i - 1] = d[position_new][i][0][0];    // Down
			}
			break;
		case AOD:
			assert(elevator.position > 1);
			position_new = elevator.position;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][0][1];  // Up
				distance[2 * i - 1] = d[position_new][i][0][0];    // Down
			}
			break;
		case AU_INV:
			// Go to floor one below and return up
			assert(elevator.elevatorState == EMPTY);
			assert(elevator.position > 1);
			position_new = elevator.position;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][1][1] + 3;  // Up
				distance[2 * i - 1] = d[position_new][i][1][0] + 3;    // Down
			}
			distance[2 * (position_new - 1)] = 1;
			break;
		case AD_INV:
			// Go to floor one above and return up
			assert(elevator.elevatorState == EMPTY);
			assert(elevator.position < N);
			position_new = elevator.position;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][0][1] + 3;  // Up
				distance[2 * i - 1] = d[position_new][i][0][0] + 3;    // Down
			}
			distance[2 * (position_new + 1) - 1] = 1;
			break;
		case AU_GR:
			assert(elevator.elevatorState == EMPTY);
			assert(elevator.position > 1);
			position_new = elevator.position;
			int groundDistance = elevator.position - 1;
			for (int i = 1; i <= N; ++i) {
				distance[2 * i] = d[position_new][i][1][1] + 2 * groundDistance + 1;  // Up
				distance[2 * i - 1] = d[position_new][i][1][0] + 2 * groundDistance + 1;    // Down
			}
			for (int j = 1; j < position_new; ++j) {
				distance[2*j] = groundDistance + 1 + j - 1;
			}
			distance[2] = groundDistance;
			break;
//		case AS:
//			// TODO
//			break;
	}
	distance[1] = 0;    // Exclude going down at 1
	distance[2 * N] = 0;    // Exclude going up at N

	return distance;
}

// cost of people inside lift after taking action
double State::insideCost(Elevator &elevator, ElevatorAction elevatorAction) {
	if ((elevator.elevatorState == JUST_FULL) or
	    (elevator.elevatorState == FULL and (elevatorAction == AU or elevatorAction == AD))) {
		return (elevator.getNumberOfPeople() * WAIT_COST);
	}
	else {
		return 0;
	}
}


double State::getMinCost(int distance1[2 * N + 1], int distance2[2 * N + 1]) {
	double cost = 0;
	int d_min_up, d_min_down;
	int up_time, down_time;
	for (int i = 1; i <= N; ++i) {
		d_min_up = min(distance1[2 * i], distance2[2 * i]);
		d_min_down = min(distance1[2 * i - 1], distance2[2 * i - 1]);
		up_time = time_up[i] + d_min_up;
		down_time = time_down[i] + d_min_down;
		assert(up_time <= T);
		assert(down_time <= T);
		// wait cost
		cost += (n_exp_up[i][up_time][down_time] * WAIT_COST);
		cost += (n_exp_down[i][up_time][down_time] * WAIT_COST);
		// electricity cost
		cost += (d_min_up * ELECTRICITY_COST);
		cost += (d_min_down * ELECTRICITY_COST);
	}

	return cost;
}

// get action according to policy
Action State::getPolicyAction() {
	// Update state if full
	if (elevator1.elevatorState == FULL)
		elevator1.updateFullState();
	if (elevator2.elevatorState == FULL)
		elevator2.updateFullState();

	// Hold all possible actions for elevators
	vector<ElevatorAction> actions1;
	vector<ElevatorAction> actions2;

	// Get actions for elevators
	actions1 = elevator1.getActions(actionq1);
	actions2 = elevator2.getActions(actionq2);

	double minCost = INT_MAX;
	Action greedyAction;
	for (auto action1: actions1) {
		for (auto action2: actions2) {
			int *distance1 = getDistanceArr(elevator1, action1);
			int *distance2 = getDistanceArr(elevator2, action2);
			// Greedy policy
			double cost = getMinCost(distance1, distance2) + insideCost(elevator1, action1) +
			              insideCost(elevator2, action2);
			if (cost < minCost) {
				greedyAction = make_pair(action1, action2);
				minCost = cost;
			}
		}
	}

	// Modify action for empty elevator
	if(elevator1.elevatorState == EMPTY) {
		if(greedyAction.first == AU_INV or greedyAction.first == AU_GR or greedyAction.first == AD_INV) {
			greedyAction.first = actionq1.front();
			actionq1.pop();
		}
	}
	if(elevator2.elevatorState == EMPTY) {
		if(greedyAction.second == AU_INV or greedyAction.second == AU_GR or greedyAction.second == AD_INV) {
			greedyAction.second = actionq2.front();
			actionq2.pop();
		}
	}

	// Update elevator states on basis of greedy action found
	// This step is not done in simulator as simulator is out of our control in real problem
	if (elevator1.elevatorState != FULL)
		elevator1.updateState(greedyAction.first, actionq1);
	if (elevator2.elevatorState != FULL)
		elevator2.updateState(greedyAction.second, actionq2);

	return greedyAction;
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
	// people inside elevator
	cost += (this->elevator1.getNumberOfPeople() * WAIT_COST);
	cost += (this->elevator2.getNumberOfPeople() * WAIT_COST);
	// people outside elevator
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= N; ++j) {
			cost += num_out[i][j];
		}
	}
	return cost;
}


double State::runSimulation(int epochs) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));
	// clear action queues
	actionq1 = queue<ElevatorAction>();
	actionq2 = queue<ElevatorAction>();

	State *startState = new State(*this);   // copy state for simulation
	// Randomly initialize left position

	Action action = startState->getPolicyAction();
	startState->simulateStep(num_out);  // simulator step
	cost += startState->applyAction(action, num_out);  // change lift state according to action, add action cost
	cost += startState->getWaitCost(action, num_out);  // Add wait cost

	return cost;
}

string State::toString() {
	string stateStr = "";
	stateStr += "Elevator1:\n";
	stateStr += (elevator1.toString() + "\n");
	stateStr += "Elevator2:\n";
	stateStr += (elevator2.toString() + "\n");
	string upTimeStr = "Up Times: ";
	for (int i = 1; i <= N; ++i) {
		upTimeStr += (to_string(time_up[i]) + " ");
	}
	stateStr += (upTimeStr + "\n");
	string downTimeStr = "Down Times: ";
	for (int i = 1; i <= N; ++i) {
		downTimeStr += (to_string(time_down[i]) + " ");
	}
	stateStr += (downTimeStr + "\n");
	return stateStr;
}

