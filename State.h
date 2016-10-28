#ifndef LIFT_AGENT_STATE_H
#define LIFT_AGENT_STATE_H

#include "utility.h"
#include <cstdlib>
#include <cstring>
#include <queue>

extern double p, q, r;
const int T = 2 * N - 1;
extern double P[T + 2][N + 1][N + 1][T + 2][T + 2];
extern double P_arrived[N + 1][T + 2][T + 2];
extern double n_exp[N + 1][N + 1][T + 2][T + 2];
extern double n_exp_up[N + 1][T + 2][T + 2];
extern double n_exp_down[N + 1][T + 2][T + 2];

extern int d[N + 1][N + 1][2][2];

class Elevator {
public:
	int position;
	int alight[N + 1];
	ElevatorState elevatorState;
	bool is_up;

	Elevator() {
		this->position = rand() % (N - 1) + 1;
		memset(alight, 0, sizeof(alight));
		elevatorState = EMPTY;
		is_up = (bool) (rand() % 2);
	}

	// check if elevator is empty
	bool isEmpty() {
		if (is_up) {
//			if(position == N)
//				return true;
			int num_alight = 0;
			for (int i = 1; i <= N; ++i) {  // TODO: Check if only position to N is required
				num_alight += alight[i];
			}
			return (num_alight == 0);
		}
		else {
//			if(position == 0)
//				return true;
			int num_alight = 0;
			for (int i = 1; i <= N; ++i) {
				num_alight += alight[i];
			}
			return (num_alight == 0);
		}
	}

	// return number of people in elevator
	int getNumberOfPeople() {
		int num = 0;
		for (int alightFloor: alight) {
			num += alightFloor;
		}
		return num;
	}

	vector<ElevatorAction> getActions(queue <ElevatorAction> &q) {
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
		else if (this->elevatorState == EMPTY) {
			if (position > 1) {
				actions.push_back(AD);
				actions.push_back(AOD);
			}
			if (position < N) {
				actions.push_back(AU);
				actions.push_back(AOU);
			}
			if(position > 1 and position < N) {
				actions.push_back(AU_INV);
				actions.push_back(AD_INV);
			}
			if(position == 3) {     // TODO: explore conditions in greater detail
				actions.push_back(AU_GR);
			}
//			actions.push_back(AS);  // TODO: add support for AS
		}
		else if (this->elevatorState == JUST_FULL) {
			assert(not q.empty());
			actions.push_back(q.front());
			q.pop();
		}
		return actions;
	}

	void updateState(ElevatorAction elevatorAction, queue <ElevatorAction> &q) {
		if (elevatorState == EMPTY) {
			updateEmptyState(elevatorAction, q);
		}
		else if (elevatorState == JUST_FULL) {
			updateJustFullState(elevatorAction, q);
		}
	}

	void updateEmptyState(ElevatorAction elevatorAction, queue <ElevatorAction> &q) {
		assert(elevatorState == EMPTY);
		if (elevatorAction == AU or elevatorAction == AOU) {
			this->elevatorState = FULL;
			is_up = true;
		}
		else if (elevatorAction == AD or elevatorAction == AOD) {
			this->elevatorState = FULL;
			is_up = false;
		}
		else if (elevatorAction == AU_INV) {
			assert(q.empty());
			q.push(AD);
			q.push(AOU);
			q.push(AU);
			this->elevatorState = JUST_FULL;
			is_up = true;
		}
		else if (elevatorAction == AD_INV) {
			assert(q.empty());
			q.push(AU);
			q.push(AOD);
			q.push(AD);
			this->elevatorState = JUST_FULL;
			is_up = false;
		}
		else if (elevatorAction == AU_GR) {
			assert(q.empty());
			for (int i = 1; i < position; ++i) {
				q.push(AD);
			}
			q.push(AOU);
			for (int i = 1; i < position; ++i) {
				q.push(AU);
			}
			this->elevatorState = JUST_FULL;
			is_up = true;
		}
	}

	void updateJustFullState(ElevatorAction elevatorAction, queue <ElevatorAction> &q) {
		assert(elevatorState == JUST_FULL);
		if (q.empty()) {
			this->elevatorState = FULL;
		}
	}

	void updateFullState() {
		assert(elevatorState == FULL);
		if (isEmpty()) {
			this->elevatorState = EMPTY;
		}
		if (is_up and position == N) {
			is_up = false;
		}
		if (not (is_up) and position == 1) {
			is_up = true;
		}

	}

	double applyAction(ElevatorAction action, int num_out[N + 1][N + 1]) {
		assert((action == AU_INV) and (action == AU_GR) and (action == AD_INV));
		double cost = 0;
		if ((action == AU) and (position < N)) {
			position++;
			cost = 1;
		}
		else if ((action == AD) and (position > 1)) {
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

	string toString() {
		string elevatorStr = "";
		elevatorStr += (to_string(position) + "\n");
		string alightStr = "Alight: ";
		for (int i = 1; i <= N; ++i) {
			alightStr += (to_string(alight[i]) + " ");
		}
		elevatorStr += (alightStr + "\n");
		string stateStr;
		switch (elevatorState) {
			case EMPTY:
				stateStr = "EMPTY";
				break;
			case JUST_FULL:
				stateStr = "JUST_FULL";
				break;
			case FULL:
				stateStr = "FULL";
				break;
		}
		elevatorStr += (stateStr + "\n");
		elevatorStr += to_string(is_up);
		return elevatorStr;
	}
};

class State {
public:
	Elevator elevator1;
	Elevator elevator2;
	int time_up[N + 1];
	int time_down[N + 1];

	State();

	State *getResState(Action action);

	int *getDistanceArr(Elevator &elevator, ElevatorAction elevatorAction);

	double insideCost(Elevator &elevator, ElevatorAction elevatorAction);

	double getMinCost(int distance1[2 * N + 2], int distance2[2 * N + 2]);

	Action getPolicyAction();

	void simulateStep(int num_out[N + 1][N + 1]);

	double applyAction(Action action, int num_out[N + 1][N + 1]);

	double getWaitCost(Action action, int num_out[N + 1][N + 1]);

	double runSimulation(int epochs);

	string toString();
};

#endif //LIFT_AGENT_STATE_H
