#include <climits>
#include <cassert>

using namespace std;

const int N = 5;
const int K = 2;

const double WAIT_COST = 2.0;
const double ELECTRICITY_COST = 1.0;

const double C_policy = 1.5;

enum ElevatorAction {
	AU,
	AD,
	AOU,
	AOD,
	AS,
	AU_INV,
	AD_INV,
	AU_GR
};

typedef pair<ElevatorAction, ElevatorAction> Action;

enum ElevatorState {
	FULL,
	EMPTY,
	JUST_FULL
};

// Combination matrix
const int C[11][11] = {
		{1, 0,  0,  0,   0,   0,   0,   0,   0,  0,  0},
		{1, 1,  0,  0,   0,   0,   0,   0,   0,  0,  0},
		{1, 2,  1,  0,   0,   0,   0,   0,   0,  0,  0},
		{1, 3,  3,  1,   0,   0,   0,   0,   0,  0,  0},
		{1, 4,  6,  4,   1,   0,   0,   0,   0,  0,  0},
		{1, 5,  10, 10,  5,   1,   0,   0,   0,  0,  0},
		{1, 6,  15, 20,  15,  6,   1,   0,   0,  0,  0},
		{1, 7,  21, 35,  35,  21,  7,   1,   0,  0,  0},
		{1, 8,  28, 56,  70,  56,  28,  8,   1,  0,  0},
		{1, 9,  36, 84,  126, 126, 84,  36,  9,  1,  0},
		{1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1}
};
