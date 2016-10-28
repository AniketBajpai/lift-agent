#ifndef LIFT_AGENT_STATE_H
#define LIFT_AGENT_STATE_H

#include "utility.h"

class Elevator {
public:
	int position;
	int alight[N];
	ElevatorState elevatorState;
	bool is_up;

	bool isTerminalPosition() {
		return ((is_up and position == N) or (not (is_up) and position == 1));
	}

	int getNumberOfPeople() {
		int num = 0;
		for (int alightFloor: alight) {
			num += alightFloor;
		}
		return num;
	}

	vector<ElevatorAction> getActions() {
		vector<ElevatorAction> actions;
		if (this->elevatorState == FULL) {
			if (is_up) {
				actions.push_back(AU);
				actions.push_back(AOU);
			}
			else if (not is_up) {
				actions.push_back(AD);
				actions.push_back(AOD);
			}
		}
		else {
			if (not (is_up and position == N)) {
				actions.push_back(AD);
				actions.push_back(AOD);
			}
			if (not (not (is_up) and position == 1)) {
				actions.push_back(AU);
				actions.push_back(AOU);
			}
			// actions.push_back(AS);
		}
		return actions;
	}

	double applyAction(ElevatorAction action, int num_out[N + 1][N + 1]) {
		double cost = 0;
		if ((action == AU) and (position <= N)) {
			position++;
			cost = 1;
		}
		else if ((action == AD) and (position >= 1)) {
			position--;
			cost = 1;
		}
		else if (action == AOU) {
			alight[position] = 0;
			for (int i = position + 1; i <= N; ++i) {
				alight[i] += num_out[position][i];
				num_out[position][i] = 0;
			}
		}
		else if (action == AOD) {
			alight[position] = 0;
			for (int i = 1; i < position; ++i) {
				alight[i] += num_out[position][i];
				num_out[position][i] = 0;
			}
		}
		return cost;
	}
};

class State {
public:
	Elevator elevator1;
	Elevator elevator2;
	int time_up[N];
	int time_down[N];

	vector<Action> getActions();

	State *getResState(Action action);

	pair<Action, State *> getRandomNextState();

	Action getPolicyAction();

	void simulateStep(int num_out[N + 1][N + 1]);

	double applyAction(Action action, int num_out[N + 1][N + 1]);

	double getWaitCost(Action action, int num_out[N + 1][N + 1]);

	double runSimulation(int epochs);
};

#endif //LIFT_AGENT_STATE_H
