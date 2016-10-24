#ifndef LIFT_AGENT_STATE_H
#define LIFT_AGENT_STATE_H

#include "utility.h"

class State
{
public:
    int position1; int position2;
    int alight1[N]; int alight2[N];
    int up_open1[N]; int up_open2[N];
    int down_open[N]; int down_open2[N];
    bitset<N> out_pressed_up;
    bitset<N> out_pressed_down;
    bool is_up;
};

#endif //LIFT_AGENT_STATE_H
