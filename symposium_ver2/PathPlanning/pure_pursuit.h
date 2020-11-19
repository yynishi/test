#ifndef _PURE_PURSUIT_H_
#define _PURE_PURSUIT_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

class State{
    public:
        float x, y, yaw, v, rear_x, rear_y;
    public:
        State();
        State(float x_, float y_, float yaw_, float v_);
        void update(float a, float delta);
        float calc_distance(float point_x, float point_y);
};


class States{
    public:
        std::vector<float> x, y, yaw, v, t;
    public:
        void append(float time, State state);
};


class TargetCourse2{
    public:
        std::vector<float> cx, cy;
        int old_nearest_point_index = -1;
    public:
        TargetCourse2(std::vector<float> cx_, std::vector<float> cy_);
        std::pair<int, float> search_target_index(State state);
};

#endif // _PURE_PURSUIT_H_