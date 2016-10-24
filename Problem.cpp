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

    for (int floor = 1; floor <= N; ++floor) {
        P_arrived[floor][0][0] = 1;     // time passed = 0
        for (int t = 1; t <= T; ++t) {  // time t
            for (int m = 0; m <= t; ++m) {
                float prob = 0;
                float arrival_choice_factor;
                if(floor == 1) {
                    arrival_choice_factor = q;
                }
                else {
                    arrival_choice_factor = (1-q)/N;
                }
                float q_factor = pow((arrival_choice_factor/(1-arrival_choice_factor)), m);
                float p_factor = p*(1-arrival_choice_factor);
                for (int k = m; k <= t; ++k) {
                    prob += (q_factor * pow(p_factor, k));
                }
                P_arrived[1][m][t] = prob;
            }
        }
    }
}

void precompute_P() {
    memset(P, 0, sizeof(P));

    for (int floor_start = 1; floor_start <= N; ++floor_start) {
        // for floor_end = floor_start, probabilities are 0
        // For t = 0, all probabilities are 0
        for (int t = 1; t <= T; ++t) {
            for (int floor_end = 2; floor_end <= N; ++floor_end) {
                for (int n = 0; n <= t; ++n) {
                    float prob = 0;
                    float n_factor, m_factor;
                    float alight_choice_factor;
                    if(floor_start == 1) {
                        alight_choice_factor = 1/(N-1);
                    }
                    else {
                        if(floor_end == 1) {
                            alight_choice_factor = r;
                        }
                        else {
                            alight_choice_factor = (1-r)/(N-2);
                        }
                    }
                    n_factor = pow(alight_choice_factor/(1-alight_choice_factor), n);
                    m_factor = 1-alight_choice_factor;
                    for (int m = n; m <= t; ++m) {
                        prob += (n_factor * pow(m_factor, m) * P_arrived[floor_start][m][t]);
                    }
                    P[n][1][floor_end][t] = prob;
                }
            }
        }
    }
}

void precompute_n_exp() {
    memset(n_exp, 0, sizeof(n_exp));

    for (int floor_start = 1; floor_start <= N; ++floor_start) {
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
    precompute_n_exp_up_down();
}

int main() {
    precompute();
}