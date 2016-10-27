#ifndef LIFT_AGENT_STATE_H
#define LIFT_AGENT_STATE_H

#include "utility.h"

class Elevator
{
public:
    int position;
    int alight[N];
    ElevatorState elevatorState;
    bool is_up;

    bool isTerminalPosition() {
        return ((is_up and position==N) or (not(is_up) and position==1));
    }
};

class State
{
public:
    Elevator elevator1;
    Elevator elevator2;
    int time_up[N];
    int time_down[N];

	vector<Action> getActions();
	State* getResState(Action action);
	pair<Action, State*> getRandomNextState();
	double runSimulation(int epochs);
};

#endif //LIFT_AGENT_STATE_H
