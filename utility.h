#include <bitset>

using namespace std;

const int N = 5;

enum Action {
    AU,
    AD,
    AOU,
    AOD,
    AS
};

enum ElevatorState {
    FULL,
    EMPTY,
    JUST_FULL
};