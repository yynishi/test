#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <matplotlibcpp.h>

namespace plt = matplotlibcpp;


float k = 0.1;  // look forward gain
float Lfc = 2.0;  // [m] look-ahead distance
float Kp = 1.0;  // speed proportional gain
float dt = 0.1;  // [s] time tick
float WB = 2.9;  // [m] wheel base of vehicle

bool show_animation = true;


class State{
    public:
        float x, y, yaw, v, rear_x, rear_y;
    public:
        State();
        State(float x_, float y_, float yaw_, float v_);
        void update(float a, float delta);
        float calc_distance(float point_x, float point_y);
};

State::State() : x(0.0) ,y(0.0), yaw(0.0), v(0.0), rear_x(x - ((WB / 2) * cos(yaw))), rear_y(y - ((WB / 2) * sin(yaw))){}
State::State(float x_, float y_, float yaw_, float v_) : x(x_), y(y_), yaw(yaw_), v(v_), rear_x(x - ((WB / 2) * cos(yaw))), rear_y(y - ((WB / 2) * sin(yaw))){}
void State::update(float a, float delta){
    x += v * cos(yaw) * dt;
    y += v * sin(yaw) * dt;
    yaw += v / WB * tan(delta) * dt;
    v += a * dt;
    rear_x = x - ((WB / 2) * cos(yaw));
    rear_y = y - ((WB / 2) * sin(yaw));
}
float State::calc_distance(float point_x, float point_y){
    float dx, dy;
    dx = rear_x - point_x;
    dy = rear_y - point_y;
    return hypot(dx, dy);
}


class States{
    public:
        std::vector<float> x, y, yaw, v, t;
    public:
        void append(float time, State state);
};

void States::append(float time, State state){
    x.push_back(state.x);
    y.push_back(state.y);
    yaw.push_back(state.yaw);
    v.push_back(state.v);
    t.push_back(time);
}


class TargetCourse2{
    public:
        std::vector<float> cx, cy;
        int old_nearest_point_index = -1;
    public:
        TargetCourse2(std::vector<float> cx_, std::vector<float> cy_);
        std::pair<int, float> search_target_index(State state);
};

TargetCourse2::TargetCourse2(std::vector<float> cx_, std::vector<float> cy_) : cx(cx_), cy(cy_){}
std::pair<int, float> TargetCourse2::search_target_index(State state){
    std::vector<float> dx, dy, d;
    float distance_this_index, distance_next_index, Lf;
    int i, ind;

    if(old_nearest_point_index == -1){
        for(auto& icx : cx){
            dx.push_back(state.rear_x - icx);
        }
        for(auto& icy : cy){
            dy.push_back(state.rear_y - icy);
        }
        for(i = 0; i < dx.size(); i++){
            d.push_back(hypot(dx[i], dy[i]));
        }
        ind = std::min_element(d.begin(),d.end()) - d.begin();
        old_nearest_point_index = ind;
    }
    else{
        ind = old_nearest_point_index;
        distance_this_index = state.calc_distance(cx[ind], cy[ind]);
        while(true){
            distance_next_index = state.calc_distance(cx[ind + 1], cy[ind + 1]);
            if(distance_this_index < distance_next_index) break;
            if((ind + 1) < cx.size()) ind++;
            distance_this_index = distance_next_index;
        }
        old_nearest_point_index = ind;
    }

    Lf = k * state.v + Lfc;  // update look ahead distance

    while(Lf > state.calc_distance(cx[ind], cy[ind])){
        if((ind + 1) >= cx.size()) break;
        ind++;
    }

    return std::make_pair(ind, Lf);
}


float proportional_control(float target, float current){
    float a;
    a = Kp * (target - current);
    return a;
}


std::pair<float, int> pure_pursuit_steer_control(State state, TargetCourse2& trajectory, int pind){
    std::pair<int, float> ind_Lf = trajectory.search_target_index(state);
    int ind = ind_Lf.first;
    float Lf = ind_Lf.second, tx, ty, alpha, delta;

    if(pind >= ind){
        ind = pind;
    }
    if(ind < trajectory.cx.size()){
        tx = trajectory.cx[ind];
        ty = trajectory.cy[ind];
    }
    else{
        tx = trajectory.cx[trajectory.cx.size() - 1];
        ty = trajectory.cy[trajectory.cy.size() - 1];
        ind = trajectory.cx.size() - 1;
    }

    alpha = atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw;

    delta = atan2(2.0 * WB * sin(alpha) / Lf, 1.0);
    
    return std::make_pair(delta, ind);
}

int main(){
    float i, ai, di;
    std::pair<float, int> di_target_ind;
    std::vector<float> cx, cy;
    //  target course
    for(i = 0; i < 50; i = i + 0.5){
        cx.push_back(i);
    }
    for(auto& ix : cx){
        cy.push_back(sin(ix / 5.0) * ix / 2.0);
    }

    float target_speed = 10.0 / 3.6;  // [m/s]

    float T = 100.0;  // max simulation time

    // initial state
    State state(0.0, -3.0, 0.0, 0.0);

    int lastIndex = cx.size() - 1;
    float time = 0.0;
    States states;
    states.append(time, state);
    TargetCourse2 target_course(cx, cy);
    int target_ind = target_course.search_target_index(state).first;

    while(T >= time & lastIndex > target_ind){
        // Calc control input
        ai = proportional_control(target_speed, state.v);
        di_target_ind = pure_pursuit_steer_control(state, target_course, target_ind);
        di = di_target_ind.first;
        target_ind = di_target_ind.second;

        // std::cout << std::fixed << std::setprecision(2) << time << " " << ai << " " << di << std::endl;
        state.update(ai, di);  // Control vehicle
        // std::cout << std::fixed << std::setprecision(2) << time << " " << state.x << " " << state.y << " " << state.yaw << " " << state.v << std::endl;

        time += dt;
        states.append(time, state);

        if(show_animation){  // pragma: no cover
            plt::cla();
            // for stopping simulation with the esc key.
            plt::plot(cx, cy, "-r");
            plt::plot(states.x, states.y, "-b");
            // plt::plot(cx[target_ind], cy[target_ind], "xg");
            plt::axis("equal");
            plt::grid(true);
            plt::title("Speed[km/h]:" + std::to_string(state.v * 3.6));
            plt::pause(0.001);
        }
    }
}