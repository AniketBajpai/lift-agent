#include <bits/stdc++.h>
#include <random>
#include "State.h"


State::State() {
	elevator1 = Elevator();
	elevator2 = Elevator();
	memset(time_up, 0, sizeof(time_up));
	memset(time_up, 0, sizeof(time_up));
	memset(time_down, 0, sizeof(time_down));
	actionq1 = queue<ElevatorAction>();
	actionq2 = queue<ElevatorAction>();
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
	int *distance = new int[2 * N + 1];
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
				distance[2 * j] = groundDistance + 1 + j - 1;
			}
			distance[2] = groundDistance;
			break;
//		case AS:
//			position_new = elevator.position;
//			for (int i = 1; i <= N; ++i) {
//				distance[2 * i] = d[position_new][i][1][1];  // Up
//				distance[2 * i - 1] = d[position_new][i][1][0];    // Down
//			}
//			break;
	}
	distance[1] = 0;    // Exclude going down at 1
	distance[2 * N] = 0;    // Exclude going up at N

	return distance;
}


int *State::getMinDistanceArr(int *distance1, int *distance2) {
	int *distance = new int[2 * N + 1];
	for (int i = 1; i <= N; ++i) {
		distance[2 * i] = min(distance1[2 * i], distance2[2 * i]);
		distance[2 * i - 1] = min(distance1[2 * i - 1], distance2[2 * i - 1]);
	}
	return distance;
}

// cost of people inside lift after taking action
double State::insideCost(Elevator &elevator, int *distance) {
	if (elevator.elevatorState == FULL) {
		double cost = 0;
		for (int i = 1; i <= N; ++i) {
			if (elevator.is_up) {
				cost += (elevator.alight[i] * distance[2 * i] * WAIT_COST);
			}
			else {
				cost += (elevator.alight[i] * distance[2 * i - 1] * WAIT_COST);
			}
		}
		return cost;
	}
	else {
		return 0;
	}
}


double State::getMinCost(int distance1[2 * N + 1], int distance2[2 * N + 1]) {
	double cost = 0;
	const double MAX_COST = N * T;
	int d_min_up, d_min_down;
	int up_time, down_time;
//	cout << "Distances: " << endl;
	for (int i = 1; i <= N; ++i) {
		d_min_up = min(distance1[2 * i], distance2[2 * i]);
		d_min_down = min(distance1[2 * i - 1], distance2[2 * i - 1]);
//		cout << d_min_up << " " << d_min_down << endl;
		if (time_up[i] != 0) {
			up_time = time_up[i] + d_min_up;
		}
		else {
			double time = 0;
			for (int j = 1; j <= T; ++j) {
				if (j > d_min_up)
					break;
				time += (P_arrived[i][i][j] * (d_min_up - j));
			}
			up_time = (int) round(time);
		}
		if (time_down[i] != 0) {
			down_time = time_down[i] + d_min_down;
		}
		else {
			double time = 0;
			for (int j = 1; j <= T; ++j) {
				if (j > d_min_down)
					break;
				time += (P_arrived[i][i][j] * (d_min_down - j));
			}
			down_time = (int) round(time);
		}
		if (up_time <= T and down_time <= T) {
			assert(up_time <= T);
			assert(down_time <= T);
			// wait cost
			cost += (n_exp_up[i][up_time][down_time] * WAIT_COST);
			cost += (n_exp_down[i][up_time][down_time] * WAIT_COST);
			// Electricity cost not counted for now, as it is inevitable
//			// electricity cost
//			cost += (d_min_up * ELECTRICITY_COST);
//			cost += (d_min_down * ELECTRICITY_COST);
		}
		else {
			cost = MAX_COST; // TODO: improve
		}
	}

	return cost;
}

// get action according to policy
pair<Action, double> State::getPolicyAction() {
	// // Update state if full
	// if (elevator1.elevatorState == FULL)
	// 	elevator1.updateFullState();
	// if (elevator2.elevatorState == FULL)
	// 	elevator2.updateFullState();

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
			double cost = getMinCost(distance1, distance2) +
			              insideCost(elevator1, distance1) +
			              insideCost(elevator2, distance2);
			if (cost < minCost) {
				greedyAction = make_pair(action1, action2);
				minCost = cost;
			}
		}
	}

	// TODO: delete action vectors and free memory

	// Modify action for empty elevator
	if (elevator1.elevatorState == EMPTY) {
		if (greedyAction.first == AU_INV or greedyAction.first == AU_GR or greedyAction.first == AD_INV) {
			greedyAction.first = actionq1.front();
			actionq1.pop();
		}
	}
	if (elevator2.elevatorState == EMPTY) {
		if (greedyAction.second == AU_INV or greedyAction.second == AU_GR or greedyAction.second == AD_INV) {
			greedyAction.second = actionq2.front();
			actionq2.pop();
		}
	}

	// // Update elevator states on basis of greedy action found
	// // This step is not done in simulator as simulator is out of our control in real problem
	// if (elevator1.elevatorState != FULL)
	// 	elevator1.updateState(greedyAction.first, actionq1);
	// if (elevator2.elevatorState != FULL)
	// 	elevator2.updateState(greedyAction.second, actionq2);

	return make_pair(greedyAction, minCost);
}

vector<pair<Action, double> > State::getActionCosts() {
	// Update state if full

	// Hold all possible actions for elevators
	vector<ElevatorAction> actions1;
	vector<ElevatorAction> actions2;

	// Get actions for elevators
	actions1 = elevator1.getActions(actionq1);
	actions2 = elevator2.getActions(actionq2);

	vector<pair<Action, double> > costs;
	for (auto action1: actions1) {
		for (auto action2: actions2) {
			int *distance1 = getDistanceArr(elevator1, action1);
			int *distance2 = getDistanceArr(elevator2, action2);
			// Greedy policy
			double cost = getMinCost(distance1, distance2) +
			              insideCost(elevator1, distance1) +
			              insideCost(elevator2, distance2);
			Action action = make_pair(action1, action2);
			costs.push_back(make_pair(action, cost));
		}
	}

	return costs;
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
		num_out[i][j]++;
	}

}

double State::applyAction(Action action, int num_out[N + 1][N + 1]) {
	double cost = 0;
	cost += this->elevator1.applyAction(action.first, num_out);
	cost += this->elevator1.applyAction(action.second, num_out);
	if(action.first == AOU) {
		time_up[elevator1.position] = 0;
	}
	else if(action.first == AOD) {
		time_down[elevator1.position] = 0;
	}
	if(action.second == AOU) {
		time_up[elevator2.position] = 0;
	}
	else if(action.second == AOD) {
		time_down[elevator2.position] = 0;
	}

	return cost;
}


double State::getWaitCost(int num_out[N + 1][N + 1]) {
	double cost = 0;
	// people inside elevator
	cost += (this->elevator1.getNumberOfPeople() * WAIT_COST);
	cost += (this->elevator2.getNumberOfPeople() * WAIT_COST);
	// people outside elevator
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= N; ++j) {
			cost += (num_out[i][j] * WAIT_COST);
		}
	}
	return cost;
}

vector<Action> State::getSimpleActions() {
	vector<Action> actions;
	vector<ElevatorAction> actions1 = elevator1.getSimpleActions();
	vector<ElevatorAction> actions2 = elevator2.getSimpleActions();
	for (auto action1: actions1) {
		for (auto action2: actions2) {
			actions.push_back(make_pair(action1, action2));
		}
	}
	return actions;
}

double runSimulation(State *startState, Action action, int epochs) {
	double cost = 0;
	int num_out[N + 1][N + 1];  // number of people going from floor1 to floor2
	memset(num_out, 0, sizeof(num_out));
//	// clear action queues
//	actionq1 = queue<ElevatorAction>();
//	actionq2 = queue<ElevatorAction>();

//	State *startState = new State(*this);   // copy state for simulation

	startState->simulateStep(num_out);
	cost += startState->applyAction(action, num_out);   // Initial action cost
	cost += startState->getWaitCost(num_out);

	ElevatorState prevState, currState;
	while (epochs > 0) {
		prevState = startState->elevator1.elevatorState;
		Action chosenAction = startState->getPolicyAction().first;
		startState->simulateStep(num_out);  // simulator step
		currState = startState->elevator1.elevatorState;
		cost += startState->applyAction(chosenAction,
		                                num_out);  // change lift state according to action, add action cost
		cost += startState->getWaitCost(num_out);  // Add wait cost
		bool isNewEpoch = (prevState == FULL and currState == EMPTY);
		if (isNewEpoch)
			epochs--;
	}

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
//helper fxns for splitting
/*void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}*/
//end of helper function for splitting 
void State::update(string s) {
	vector<string> input; //= split("s", ' ');	

	istringstream iss(s);
	copy(istream_iterator<string>(iss),
	     istream_iterator<string>(),
	     back_inserter(input));

	//cout<<input[0]<<endl;
	//cout<<input.size(); 
	if (input[0] == "0" && input.size() == 1)
		return;
	for (int i = 0; i < input.size(); i++) {
		char a = input[i][1];
		char b = input[i][3];
		int f = b - 48;
		//cout<<a<<b;
		if (a == 'U') {
			time_up[f]++;
		}
		else if (a == 'D') {
			time_down[f]++;
		}
		else {
			if (input[i][4] == '1') {    //cout<<"yo";
				elevator1.btnPressed[input[i][2] -48] = true;    //TODO: make button pressed false when liy opens at a floor
				elevator1.alight[input[i][2] -48] = 5;
			}
			else if (input[i][4] == '2') {    //cout<<"yo2";
				elevator2.btnPressed[input[i][2] - 48] = true;
				elevator2.alight[input[i][2] -48] = 5;
			}
		}

	}
}

void printS(State s) {
	cerr<<"elev1: "<<s.elevator1.position<<" elev2: "<<s.elevator2.position<<endl;
	cerr << "TIME UP:";
	for (int i = 1; i < N + 1; i++) {
		cerr << s.time_up[i] << " ";
	}
	cerr << endl;
	cerr << "TIME DOWN:";
	for (int i = 1; i < N + 1; i++) {
		cerr << s.time_down[i] << " ";
	}
	cerr << endl;

	cerr << "ALIGHT elevator1" << endl;
	for (int i = 1; i < N + 1; i++) {
		cerr << s.elevator1.alight[i] << " ";
	}
	cerr << endl;

	cerr << "BTN pressed elevator1" << endl;
	for (int i = 1; i < N + 1; i++) {
		cerr << s.elevator1.btnPressed[i] << " ";
	}
	cerr << endl;

	cerr << "ALIGHT elevator2" << endl;
	for (int i = 1; i < N + 1; i++) {
		cerr << s.elevator2.alight[i] << " ";
	}
	cerr << endl;
	cerr << "BTN pressed elevator2" << endl;
	for (int i = 1; i < N + 1; i++) {
		cerr << s.elevator2.btnPressed[i] << " ";
	}
	cerr << endl;
} 