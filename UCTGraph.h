#ifndef LIFT_AGENT_UCTGRAPH_H
#define LIFT_AGENT_UCTGRAPH_H

#include "UCTGraphNode.h"

class UCTGraph {
public:
	static vector<UCTGraphNode *> nodeList;

	static UCTGraphNode *getNodeFromState(State *state);

	static void run(State *startState, int depth);

	static double runhelper(State *state, int depth, int num_out[N + 1][N + 1]);

	static double runSimulation(State *state, Action action, int epochs);

	static vector<pair<Action, double> > getBaseCosts(State *state, int depth);
};

#endif //LIFT_AGENT_UCTGRAPH_H
