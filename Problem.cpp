#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include "State.h"

double p, q, r;

double P[
		T + 2][
		N + 1][
		N + 1][
		T + 2][T + 2];  // prob. of n people going from f1 to f2 after t1 - up waiting time, t2 - down waiting time
double P_arrived[N + 1][T + 2][T + 2];
double n_exp[
		N + 1][
		N + 1][
		T + 2][
		T + 2];   // expected no. of people going from f1 to f2 after t1 - up waiting time, t2 - down waiting time
double n_exp_up[N + 1][T + 2][T + 2];     // expected no. of people going up from f
double n_exp_down[N + 1][T + 2][T + 2];   // expected no. of people going down from f

int d[N + 1][N + 1][2][2];  // distance between elevator and person given direction of elevator and person

void precompute_P_arrival(double p, double q, double r) {
	memset(P_arrived, 0, sizeof(P_arrived));

	for (int floor = 1; floor <= N; ++floor) {
		P_arrived[floor][0][0] = 1;     // time passed = 0
		for (int t = 1; t <= T; ++t) {  // time t
			for (int m = 0; m <= t; ++m) {
				double prob = 0;
				double arrival_choice_factor;
				if (floor == 1) {
					arrival_choice_factor = q;
				}
				else {
					arrival_choice_factor = (1.0 - q) / (double) N;
				}
				for (int k = m; k <= t; ++k) {
					prob += ((double) C[t][k] * pow(p, k) * pow((1.0 - p), (t - k)) * (double) C[k][m] *
					         pow(arrival_choice_factor, m) * pow((1.0 - arrival_choice_factor), (k - m)));
				}
				P_arrived[floor][m][t] = prob;
			}
		}
	}
}

void precompute_P(double p, double q, double r) {
	memset(P, 0, sizeof(P));
	for (int t1 = 0; t1 <= T; ++t1) {
		for (int t2 = 0; t2 <= T; ++t2) {
			for (int floor_start = 1; floor_start <= N; ++floor_start) {
				// for floor_end = floor_start, probabilities are 0
				// For t = 0, all probabilities are 0
				for (int floor_end = 1; floor_end <= N; ++floor_end) {
					if (floor_end == floor_start) {
						continue;
					}

					double prob;
					double alight_choice_factor;
					if (floor_start == 1) {
						if (t1 != 0) {
							alight_choice_factor = 1.0 / (double) (N - 1);
						}
						else {
							alight_choice_factor = 0.0;
						}
					}
					else {
						if (t1 == 0 and t2 == 0) {
							alight_choice_factor = 0.0;
						}
						else if (t1 == 0 and t2 != 0) {
							if (floor_end > floor_start) {
								alight_choice_factor = 0.0;
							}
							else if (floor_end < floor_start) {
								if (floor_end == 1) {
									alight_choice_factor =
											r / (r + (double) (floor_start - 1) * ((1.0 - r) / (N - 2.0)));
								}
								else {
									alight_choice_factor = ((1.0 - r) / (N - 2.0)) /
									                       (r + (double) (floor_start - 1) * ((1.0 - r) / (N - 2.0)));
								}
							}
						}
						else if (t1 != 0 and t2 == 0) {
							if (floor_end > floor_start) {
								alight_choice_factor = 1.0 / (double) (N - floor_end);
							}
							else if (floor_end < floor_start) {
								alight_choice_factor = 0.0;
							}
						}
						else if (t1 != 0 and t2 != 0) {
							if (floor_end == 1) {
								alight_choice_factor = r;
							}
							else {
								alight_choice_factor = (1.0 - r) / (N - 2.0);
							}
						}
					}

					if (floor_end > floor_start) {   // Up
						for (int n = 0; n <= t1; ++n) {
							prob = 0;
							for (int m = max(1, n); m <= (t1 - 1) + 1; ++m) {
								prob += ((double) C[m][n] * pow(alight_choice_factor, n) *
								         pow(1.0 - alight_choice_factor, (m - n)) *
								         P_arrived[floor_start][m - 1][max(t1 - 1, 0)]);
							}
							P[n][floor_start][floor_end][t1][t2] = prob;
						}
					}
					else if (floor_end < floor_start) {  // Down
						for (int n = 0; n <= t2; ++n) {
							prob = 0;
							for (int m = max(1, n); m <= (t2 - 1) + 1; ++m) {
								prob += ((double) C[m][n] * pow(alight_choice_factor, n) *
								         pow(1.0 - alight_choice_factor, (m - n)) *
								         P_arrived[floor_start][m - 1][max(t2 - 1, 0)]);
							}
							P[n][floor_start][floor_end][t1][t2] = prob;
						}
					}
				}
			}
		}
	}
}

void precompute_n_exp(double p, double q, double r) {
	memset(n_exp, 0, sizeof(n_exp));

	for (int t1 = 0; t1 <= T; ++t1) {
		for (int t2 = 0; t2 <= T; ++t2) {
			for (int floor_start = 1; floor_start <= N; ++floor_start) {
				// for floor_end = floor_start, n_exp = 0
				for (int floor_end = 1; floor_end <= N; ++floor_end) {
					if (floor_start == floor_end) {
						continue;
					}
					if (floor_end > floor_start) {
						double exp = 0;
						for (int n = 1; n <= (t1 - 1) + 1; ++n) {
							exp += ((double) (n) * P[n][floor_start][floor_end][t1][t2]);
						}
						n_exp[floor_start][floor_end][t1][t2] = exp;
					}
					else if (floor_end < floor_start) {
						double exp = 0;
						for (int n = 1; n <= (t2 - 1) + 1; ++n) {
							exp += ((double) (n) * P[n][floor_start][floor_end][t1][t2]);
						}
						n_exp[floor_start][floor_end][t1][t2] = exp;
					}
				}
			}
		}
	}
}

void precompute_n_exp_up_down() {
	memset(n_exp_up, 0, sizeof(n_exp_up));
	memset(n_exp_down, 0, sizeof(n_exp_down));

	for (int floor_start = 1; floor_start <= N; ++floor_start) {
		for (int t1 = 0; t1 <= T; ++t1) {
			for (int t2 = 0; t2 <= T; ++t2) {
				// Down computation
				double exp_down = 0;
				for (int floor_end = 1; floor_end < floor_start; ++floor_end) {
					exp_down += n_exp[floor_start][floor_end][t1][t2];
				}
				n_exp_down[floor_start][t1][t2] = exp_down;

				// Up computation
				double exp_up = 0;
				for (int floor_end = floor_start + 1; floor_end <= N; ++floor_end) {
					exp_up += n_exp[floor_start][floor_end][t1][t2];
				}
				n_exp_up[floor_start][t1][t2] = exp_up;
			}
		}
	}
}

void precompute_distance() {
	memset(d, 0, sizeof(d));
	for (int i = 1; i <= N; ++i) {  // elevator position
		// elevator up, person up
		for (int j = 1; j <= i - 1; ++j) {
			d[i][j][1][1] = N - i + N + j;
		}
		for (int j = i; j <= N; ++j) {
			d[i][j][1][1] = (j - i);
		}

		// elevator up, person down
		for (int j = 1; j <= N; ++j) {
			d[i][j][1][0] = N - i + N - j;
		}

		// elevator down, person up
		for (int j = 1; j <= N; ++j) {
			d[i][j][0][1] = i + j;
		}

		// elevator down, person down
		for (int j = 1; j <= i; ++j) {
			d[i][j][0][0] = i - j;
		}
		for (int j = i + 1; j <= N; ++j) {
			d[i][j][0][0] = i + N + N - j;
		}
	}
}

void precompute(double p, double q, double r) {
	precompute_P_arrival(p, q, r);
	precompute_P(p, q, r);
	precompute_n_exp(p, q, r);
	precompute_n_exp_up_down();
	precompute_distance();
}

double getRemCost(Elevator elevator) {
	double rem_cost = 0.0;
	if (elevator.is_up and elevator.elevatorState != EMPTY) {
		for (int i = elevator.position + 1; i <= N; ++i) {
			rem_cost += ((double) (i - elevator.position) * elevator.alight[i]);
		}
	}
	else if (!elevator.is_up and elevator.elevatorState != EMPTY) {
		for (int i = elevator.position - 1; i >= 1; --i) {
			rem_cost += (double) ((elevator.position - i) * elevator.alight[i]);
		}
	}
	return rem_cost;
}

double evaluateTerminalState(State state) {
	// cost of remaining people in elevator
	double rem_cost = 0;
	rem_cost += getRemCost(state.elevator1);
	rem_cost += getRemCost(state.elevator2);

	// cost of waiting people
	double wait_cost = 0;
	for (int i = 1; i <= N; ++i) {

	}

	return (rem_cost + wait_cost);
}

//// returns cost of performing given action sequences on elevators in current state
//double getSequenceCost(State state, vector<ElevatorAction> action_seq1, vector<ElevatorAction> action_seq2) {
//	assert(action_seq1.size() == action_seq2.size());   // action sequences of same size
//	int num_actions = action_seq1.size();
//	for (int i = 0; i < num_actions; ++i) {
//
//	}
//}

int main() {
	p = 0.8;
	q = 0.5;
	r = 0.5;
	precompute(p, q, r);

	cout << n_exp_up[3][1][0] << endl;
	cout << n_exp_down[3][1][0] << endl;

}