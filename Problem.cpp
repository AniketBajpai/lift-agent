#include <iostream>
#include <cmath>
#include <cstring>
#include "utility.h"

const int T = 2*N-1;
double P[T+2][N+1][N+1][T+2];
double P_arrived[N+1][T+2][T+2];
double n_exp[N+1][N+1][T+2];
double n_exp_up[N+1][T+2];
double n_exp_down[N+1][T+2];

double p, q, r;
ElevatorState elevatorState;

void precompute_P_arrival() {
    memset(P_arrived, 0, sizeof(P_arrived));

    for (int floor = 1; floor <= N; ++floor) {
        P_arrived[floor][0][0] = 1;     // time passed = 0
        for (int t = 1; t <= T; ++t) {  // time t
            for (int m = 0; m <= t; ++m) {
                double prob = 0;
                double arrival_choice_factor;
                if(floor == 1) {
                    arrival_choice_factor = q;
                }
                else {
                    arrival_choice_factor = (1.0-q)/(double)N;
                }
                for (int k = m; k <= t; ++k) {
                    prob += ((double)C[t][k] * pow(p, k) * pow((1.0-p), (t-k)) * (double)C[k][m] * pow(arrival_choice_factor, m) * pow((1.0-arrival_choice_factor), (k-m)));
                }
                P_arrived[floor][m][t] = prob;
            }
        }
    }
}

void precompute_P() {
    memset(P, 0, sizeof(P));
    for (int t = 1; t <= T; ++t) {
        for (int floor_start = 1; floor_start <= N; ++floor_start) {
            // for floor_end = floor_start, probabilities are 0
            // For t = 0, all probabilities are 0
            for (int floor_end = 1; floor_end <= N; ++floor_end) {
                if(floor_end == floor_start) {
                    continue;
                }
                for (int n = 0; n <= t; ++n) {
                    double prob = 0;
                    double alight_choice_factor;
                    if(floor_start == 1) {
                        alight_choice_factor = 1.0/(double)(N-1);
                    }
                    else {
                        if(floor_end == 1) {
                            alight_choice_factor = r;
                        }
                        else {
                            alight_choice_factor = (1.0-r)/(double)(N-2);
                        }
                    }
                    for (int m = max(1,n); m <= t+1; ++m) {
                        prob += ((double)C[m][n] * pow(alight_choice_factor, n) * pow(1.0-alight_choice_factor, (m-n)) * P_arrived[floor_start][m-1][t]);
                    }
                    P[n][floor_start][floor_end][t] = prob;
                }
            }
        }
    }
}

void precompute_n_exp() {
    memset(n_exp, 0, sizeof(n_exp));
    // t = 0
    // 1st floor
    for (int floor_end = 2; floor_end <= N; ++floor_end) {
        n_exp[1][floor_end][0] = 1.0/(double)(N-1);
    }
    // Other floors
    for (int floor_start = 2; floor_start <= N; ++floor_start) {
        for (int floor_end = 1; floor_end <= N; ++floor_end) {
            if(floor_end == floor_start) {
                continue;
            }
            if(floor_end == 1) {
                n_exp[floor_start][floor_end][0] = r;
            }
            else {
                n_exp[floor_start][floor_end][0] = (1.0-r)/(double)(N-2);
            }
        }
    }

    // Later t
    for (int t = 1; t <= T; ++t) {
        for (int floor_start = 1; floor_start <= N; ++floor_start) {
            // for floor_end = floor_start, n_exp = 0
            for (int floor_end = 1; floor_end <= N; ++floor_end) {
                if(floor_start == floor_end) {
                    continue;
                }
                double exp = 0;
                for (int n = 1; n <= t+1; ++n) {
                    exp += ((double)(n) * P[n][floor_start][floor_end][t]);
                }
                n_exp[floor_start][floor_end][t] = exp;
            }
        }
    }
}

void precompute_n_exp_up_down() {
    memset(n_exp_up, 0, sizeof(n_exp_up));
    memset(n_exp_down, 0, sizeof(n_exp_down));

    for (int floor_start = 1; floor_start <= N; ++floor_start) {
        for (int t = 0; t <= T; ++t) {
            // Down computation
            double exp_down = 0;
            for (int floor_end = 1; floor_end < floor_start; ++floor_end) {
                exp_down += n_exp[floor_start][floor_end][t];
            }
            n_exp_down[floor_start][t] = exp_down;

            // Up computation
            double exp_up = 0;
            for (int floor_end = floor_start+1; floor_end <= N; ++floor_end) {
                exp_up += n_exp[floor_start][floor_end][t];
            }
            n_exp_up[floor_start][t] = exp_up;
        }
    }
}

void precompute() {
    precompute_P_arrival();
    precompute_P();
    precompute_n_exp();
    precompute_n_exp_up_down();
}

void print_exp() {
    for (int t = 0; t <= T; ++t) {
        cout << "TIME: " << t << endl;
        for (int n = 1; n <= N; ++n) {
            cout << "FLOOR: " << n << endl;
            cout << "Expected Up: " << n_exp_up[n][t] << endl;
            cout << "Expected Down: " << n_exp_down[n][t] << endl;
        }
        cout << endl << endl;
    }
}

int main() {
    p = 0.8;
    q = 0.5;
    r = 0.5;
    precompute();

//    print_exp();

}