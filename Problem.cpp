#include <iostream>
#include <cmath>
#include <cstring>
#include "utility.h"

const int T = 2*N;
float P[T][N][N][T];
float P_arrived[N][T+1][T+1];
float n_exp[N][N][T];
float n_exp_up[N][T];
float n_exp_down[N][T];

float p, q, r;
ElevatorState elevatorState;

void precompute_P_arrival() {
    memset(P_arrived, 0, sizeof(P_arrived));
    int floor;
    // 1st floor
    floor = 1;
    P_arrived[1][0][0] = 1;     // time passed = 0
    // time t
    for (int t = 1; t <= T; ++t) {
        for (int m = 0; m <= t; ++m) {
            float prob = 0;
            float q_factor = pow((q/(1-q)), m);
            float p_factor = p*(1-q);
            for (int k = m; k <= t; ++k) {
                prob += (q_factor * pow(p_factor, k));
            }
            P_arrived[1][m][t] = prob;
        }
    }
    // Other floors
    for (floor = 2; floor <= N; ++floor) {

    }
}

void precompute_P() {
    memset(P, 0, sizeof(P));
    int floor_start, floor_end;
    // 1st floor
    floor_start = 1;
    // for floor_end = floor_start, probabilities are 0
    // For t = 0, all probabilities are 0
    for (int t = 1; t <= T; ++t) {
        for (floor_end = 2; floor_end <= N; ++floor_end) {
            for (int n = 0; n <= t; ++n) {
                float prob = 0;
                float n_factor = pow((1/(N-2)), n);
                for (int m = n; m <= t; ++m) {
                    prob += (n_factor * pow((N-2)/(N-1), m) * P_arrived[1][m][t]);
                }
                P[n][1][floor_end][t] = prob;
            }
        }
    }

    // Other floors
}

void precompute_n_exp() {
    memset(n_exp, 0, sizeof(n_exp));
    int floor_start;
    // 1st floor
    floor_start = 1;
    // for floor_end = floor_start, n_exp = 0
    for (int t = 1; t <= T; ++t) {
        for (int floor_end = 2; floor_end <= N; ++floor_end) {
            float exp = 0;
            for (int n = 0; n <= t; ++n) {
                exp += ((n+1) * P[n][1][floor_end][t]);
            }
            n_exp[1][floor_end][t] = exp;
        }
    }
}

void precompute_n_exp_up_down() {
    memset(n_exp_up, 0, sizeof(n_exp_up));
    memset(n_exp_down, 0, sizeof(n_exp_down));
    
    for (int floor_start = 1; floor_start <= N; ++floor_start) {
        for (int t = 1; t <= T; ++t) {
            // Up computation
            float exp_up = 0;
            for (int floor_end = 1; floor_end < floor_start; ++floor_end) {
                exp_up += n_exp[floor_start][floor_end][t];
            }
            n_exp_up[floor_start][t] = exp_up;

            // Down computation
            float exp_down = 0;
            for (int floor_end = floor_start+1; floor_end <= N; ++floor_end) {
                exp_down += n_exp[floor_start][floor_end][t];
            }
            n_exp_down[floor_start][t] = exp_down;
        }
    }
}

void precompute() {
    precompute_P_arrival();
    precompute_P();
    precompute_n_exp();
}

int main() {
    precompute();
}