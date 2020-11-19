#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <string>
#include <numeric>
#include <limits>
#include <optional>
#include <chrono>
#include <cmath>
#include <matplotlibcpp.h>

#include "quintic_polynomials_planner.h"
#include "cubic_spline_planner.h"

namespace plt = matplotlibcpp;


int SIM_LOOP = 500;

// Parameter
float MAX_SPEED = 50.0 / 3.6;  // maximum speed [m/s]
float MAX_ACCEL = 2.0;  // maximum acceleration [m/ss]
float MAX_CURVATURE = 1.0;  // maximum curvature [1/m]
// float MAX_ROAD_WIDTH = 7.0;  // maximum road width [m]
float MAX_ROAD_WIDTH = 0.7;
// float D_ROAD_W = 1.0;  // road width sampling length [m]
float D_ROAD_W = 0.1;
float DT = 0.2;  // time tick [s]
float MAX_T = 5.0;  // max prediction time [m]
float MIN_T = 4.0;  // min prediction time [m]
// float START_SPEED = 10.0 / 3.6;
// float TARGET_SPEED = 20.0 / 3.6;  // target speed [m/s]
// float D_T_S = 5.0 / 3.6;  // target speed sampling length [m/s]
float START_SPEED = 0.5;
float TARGET_SPEED = 0.8;  // target speed [m/s]
float D_T_S = 0.2;  // target speed sampling length [m/s]
int N_S_SAMPLE = 1;  // sampling number of target speed
// float ROBOT_RADIUS = 2.0;  // robot radius [m]
float ROBOT_RADIUS = 0.2;

//car size
// float FULL_LENGTH = 4.8; // [m]
// float FULL_WIDTH = 1.7; // [m]
float FULL_LENGTH = 0.431; // [m]
float FULL_WIDTH = 0.19; // [m]

//obstacle size
float OB_FULL_LENGTH = 0.2; //[m]
float OB_FULL_WIDTH = 0.2; //[m]

//pure_pursuit用
float k = 0.1;  // look forward gain
// float Lfc = 2.0;  // [m] look-ahead distance
float Kp = 1.0;  // speed proportional gain
float dt = 0.1;  // [s] time tick
// float WB = 2.9;  // [m] wheel base of vehicle
float Lfc = 0.2;
float WB = 0.125;
int UPDATE_PATH_IND = 3;
int UPDATE_STEP = 8;
int OB_POS = 20;
float OB_SPEED = 0.5;

// cost weights
float K_J = 0.1;
float K_T = 0.1;
float K_D = 1.0;
float K_LAT = 1.0;
float K_LON = 1.0;

bool show_animation = false;


class QuarticPolynomial{
    public:
        float a0;
        float a1;
        float a2;
        float a3;
        float a4;
    public:
        QuarticPolynomial(){};
        QuarticPolynomial(float xs, float vxs, float axs, float vxe, float axe, float time);
        float calc_point(float t);
        float calc_first_derivative(float t);
        float calc_second_derivative(float t);
        float calc_third_derivative(float t);
};

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
        std::vector<float> x, y, yaw, v;
    public:
        void append(State state);
        void clear();
};

class FrenetPath{
    public:
        std::vector<float> t;
        std::vector<float> d;
        std::vector<float> d_d;
        std::vector<float> d_dd;
        std::vector<float> d_ddd;
        std::vector<float> s;
        std::vector<float> s_d;
        std::vector<float> s_dd;
        std::vector<float> s_ddd;
        float cd = 0.0;
        float cv = 0.0;
        float cf = 0.0;
        std::vector<float> x;
        std::vector<float> y;
        std::vector<float> yaw;
        std::vector<float> ds;
        std::vector<float> c;
    public:
        FrenetPath(){};
        std::pair<int, float> search_target_index(State state, int old_nearest_point_index);
};

std::pair<float, int> pure_pursuit_steer_control(State state, FrenetPath frenet_path, int old_nearest_point_index){
    std::pair<int, float> ind_Lf = frenet_path.search_target_index(state, old_nearest_point_index);
    int ind = ind_Lf.first;
    float Lf = ind_Lf.second, tx, ty, alpha, delta;

    tx = frenet_path.x[ind];
    ty = frenet_path.y[ind];

    alpha = atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw;

    delta = atan2(2.0 * WB * sin(alpha) / Lf, 1.0);
    
    return std::make_pair(delta, ind);
}

class TargetCourse{
    public:
        std::vector<float> rx;
        std::vector<float> ry;
        std::vector<float> ryaw;
        std::vector<float> rk;
        float length;
        Spline2D csp;
    public:
        TargetCourse(VectorXd x, VectorXd y);
        std::pair<int, float> search_target_index(State state, int old_nearest_point_index);
};

std::pair<float, int> pure_pursuit_steer_control(State state, TargetCourse target_course, int old_nearest_point_index){
    std::pair<int, float> ind_Lf = target_course.search_target_index(state, old_nearest_point_index);
    int ind = ind_Lf.first;
    float Lf = ind_Lf.second, tx, ty, alpha, delta;

    tx = target_course.rx[ind];
    ty = target_course.ry[ind];

    alpha = atan2(ty - state.rear_y, tx - state.rear_x) - state.yaw;

    delta = atan2(2.0 * WB * sin(alpha) / Lf, 1.0);
    
    return std::make_pair(delta, ind);
}

float proportional_control(float target, float current){
    float a;
    a = Kp * (target - current);
    return a;
}


QuarticPolynomial::QuarticPolynomial(float xs, float vxs, float axs, float vxe, float axe, float time){
    a0 = xs;
    a1 = vxs;
    a2 = axs /2.f;

    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    A << 3 * pow(time,2), 4 * pow(time, 3), 6 * time, 12 * pow(time, 2);
    b << vxe - a1 - 2 * a2 * time, axe - 2 * a2;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);

    a3 = x(0);
    a4 = x(1);
}

float QuarticPolynomial::calc_point(float t){
    float xt;
    xt = a0 + a1 * t + a2 * pow(t, 2) + a3 * pow(t, 3) + a4 * pow(t, 4);
    return xt;
}

float QuarticPolynomial::calc_first_derivative(float t){
    float xt;
    xt = a1 + 2 * a2 * t + 3 * a3 * pow(t, 2) + 4 * a4 * pow(t, 3);
    return xt;
}

float QuarticPolynomial::calc_second_derivative(float t){
    float xt;
    xt = 2 * a2 + 6 * a3 * t + 12 * a4 * pow(t, 2);
    return xt;
}

float QuarticPolynomial::calc_third_derivative(float t){
    float xt;
    xt = 6 * a3 + 24 * a4 * t;
    return xt;
}

std::vector<FrenetPath> calc_frenet_paths(float c_speed, float c_d, float c_d_d, float c_d_dd, float s0){
    std::vector<FrenetPath> frenet_paths;
    FrenetPath fp, tfp;
    QuinticPolynomial lat_qp;
    QuarticPolynomial lon_qp;
    float di, Ti, tv, t, Jp, Js, ds;

    for(di = -MAX_ROAD_WIDTH; di <= 0.1; di += D_ROAD_W){
        for(Ti = MIN_T; Ti <= MAX_T; Ti += DT){
            fp = FrenetPath();
            lat_qp = QuinticPolynomial(c_d, c_d_d, c_d_dd, di, 0.0, 0.0, Ti);

            for(t = 0; t < Ti; t += DT){
                fp.t.push_back(t);
                fp.d.push_back(lat_qp.calc_point(t));
                fp.d_d.push_back(lat_qp.calc_first_derivative(t));
                fp.d_dd.push_back(lat_qp.calc_second_derivative(t));
                fp.d_ddd.push_back(lat_qp.calc_third_derivative(t));
            }

            for(tv = TARGET_SPEED - D_T_S * N_S_SAMPLE; tv <= TARGET_SPEED + D_T_S * N_S_SAMPLE; tv += D_T_S){
                tfp = FrenetPath();
                copy(fp.t.begin(), fp.t.end(), std::back_inserter(tfp.t));
                copy(fp.d.begin(), fp.d.end(), std::back_inserter(tfp.d));
                copy(fp.d_d.begin(), fp.d_d.end(), std::back_inserter(tfp.d_d));
                copy(fp.d_dd.begin(), fp.d_dd.end(), std::back_inserter(tfp.d_dd));
                copy(fp.d_ddd.begin(), fp.d_ddd.end(), std::back_inserter(tfp.d_ddd));
                lon_qp = QuarticPolynomial(s0, c_speed, 0.0, tv, 0.0, Ti);

                for(t = 0; t < Ti; t += DT){
                    tfp.s.push_back(lon_qp.calc_point(t));
                    tfp.s_d.push_back(lon_qp.calc_first_derivative(t));
                    tfp.s_dd.push_back(lon_qp.calc_second_derivative(t));
                    tfp.s_ddd.push_back(lon_qp.calc_third_derivative(t));
                }

                Jp = std::inner_product(tfp.d_ddd.begin(), tfp.d_ddd.end(), tfp.d_ddd.begin(), 0.0);
                Js = std::inner_product(tfp.s_ddd.begin(), tfp.s_ddd.end(), tfp.s_ddd.begin(), 0.0);

                ds = pow((TARGET_SPEED - tfp.s_d.back()), 2);

                tfp.cd = K_J * Jp + K_T * Ti + K_D * pow(tfp.d.back(), 2);
                tfp.cv = K_J * Js + K_T * Ti + K_D * ds;
                tfp.cf = K_LAT * tfp.cd + K_LON * tfp.cv;

                frenet_paths.push_back(tfp);
            }
        }
    }
    return frenet_paths;
}

std::vector<FrenetPath> calc_global_paths(std::vector<FrenetPath> fplist, Spline2D csp){
    int i;
    float ix, iy, i_yaw, di, fx, fy, dx, dy;
    Eigen::Vector2d vec;

    for(auto& fp : fplist){

        for(i = 0; i < fp.s.size(); i++){
            vec = csp.calc_position(fp.s[i]);
            ix = vec(0);
            iy = vec(1);
            if(std::isnan(ix)) break;
            i_yaw = csp.calc_yaw(fp.s[i]);
            di = fp.d[i];
            fx = ix + di * cos(i_yaw + M_PI /2.f);
            fy = iy + di * sin(i_yaw + M_PI /2.f);
            fp.x.push_back(fx);
            fp.y.push_back(fy);
        }

        for(i = 0; i < fp.x.size() - 1; i++){
            dx = fp.x[i + 1] - fp.x[i];
            dy = fp.y[i + 1] - fp.y[i];
            fp.yaw.push_back(atan2(dy, dx));
            fp.ds.push_back(hypot(dx, dy));
        }

        fp.yaw.push_back(fp.yaw.back());
        fp.ds.push_back(fp.ds.back());

        for(i = 0; i < fp.yaw.size() - 1; i++){
            if(fp.yaw[i + 1] < 0 && fp.yaw[i] > 0 && fp.yaw[i] > M_PI /2.f){
                fp.c.push_back((2 * M_PI + fp.yaw[i + 1] - fp.yaw[i]) / fp.ds[i]);
                continue;
            }
            else if(fp.yaw[i + 1] > 0 && fp.yaw[i] < 0 && fp.yaw[i + 1] > M_PI /2.f){
                fp.c.push_back((-2 * M_PI + fp.yaw[i + 1] - fp.yaw[i]) / fp.ds[i]);
                continue;
            }
            fp.c.push_back((fp.yaw[i + 1] - fp.yaw[i]) / fp.ds[i]);
        }
    }

    return fplist;
}

bool check_collision(FrenetPath fp, std::vector<std::vector<float>> ob, States low_speed_car_states, TargetCourse target_course, int low_speed_car_target_ind){
    int i, j;
    float distance, min_distance, di;
    std::vector<float> cross_product;
    std::vector<std::vector<float>> next_setp_ob, car_vertex, car_rotate_vertex, ob_vertex, ob_rotate_vertex, low_speed_car_vertex, low_speed_car_rotate_vertex;
    State low_speed_car_state;
    std::pair<float, int> di_target_ind;

    //矩形で衝突判定(停止障害物)
    float ob_theta = 0.0;
    for(i = 0; i < ob.size(); i++){
        for(j = 0; j < fp.x.size(); j++){
            car_vertex = {{fp.x[j] - OB_FULL_LENGTH /2.f, fp.y[j] - OB_FULL_WIDTH /2.f},
                        {fp.x[j] - OB_FULL_LENGTH /2.f, fp.y[j] + OB_FULL_WIDTH /2.f},
                        {fp.x[j] + OB_FULL_LENGTH /2.f, fp.y[j] + OB_FULL_WIDTH /2.f}, 
                        {fp.x[j] + OB_FULL_LENGTH /2.f, fp.y[j] - OB_FULL_WIDTH /2.f}};
            car_rotate_vertex = {{(car_vertex[0][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[0][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[0][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[0][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[1][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[1][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[1][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[1][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[2][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[2][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[2][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[2][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[3][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[3][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[3][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[3][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]}};
            ob_vertex = {{ob[i][0] - OB_FULL_LENGTH /2.f, ob[i][1] - OB_FULL_WIDTH /2.f},
                        {ob[i][0] - OB_FULL_LENGTH /2.f, ob[i][1] + OB_FULL_WIDTH /2.f},
                        {ob[i][0] + OB_FULL_LENGTH /2.f, ob[i][1] + OB_FULL_WIDTH /2.f}, 
                        {ob[i][0] + OB_FULL_LENGTH /2.f, ob[i][1] - OB_FULL_WIDTH /2.f}};
            ob_rotate_vertex = {{(ob_vertex[0][0] - ob[i][0]) * cos(ob_theta) - (ob_vertex[0][1] - ob[i][0]) * sin(ob_theta) + ob[i][0], (ob_vertex[0][0] - ob[i][1]) * sin(ob_theta) + (ob_vertex[0][1] - ob[i][1]) * cos(ob_theta) + ob[i][1]},
                                {(ob_vertex[1][0] - ob[i][0]) * cos(ob_theta) - (ob_vertex[1][1] - ob[i][0]) * sin(ob_theta) + ob[i][0], (ob_vertex[1][0] - ob[i][1]) * sin(ob_theta) + (ob_vertex[1][1] - ob[i][1]) * cos(ob_theta) + ob[i][1]},
                                {(ob_vertex[2][0] - ob[i][0]) * cos(ob_theta) - (ob_vertex[2][1] - ob[i][0]) * sin(ob_theta) + ob[i][0], (ob_vertex[2][0] - ob[i][1]) * sin(ob_theta) + (ob_vertex[2][1] - ob[i][1]) * cos(ob_theta) + ob[i][1]},
                                {(ob_vertex[3][0] - ob[i][0]) * cos(ob_theta) - (ob_vertex[3][1] - ob[i][0]) * sin(ob_theta) + ob[i][0], (ob_vertex[3][0] - ob[i][1]) * sin(ob_theta) + (ob_vertex[3][1] - ob[i][1]) * cos(ob_theta) + ob[i][1]}};
            for(auto& car_ver : car_rotate_vertex){
                cross_product = {{(ob_rotate_vertex[1][0] - ob_rotate_vertex[0][0]) * (car_ver[1] - ob_rotate_vertex[0][1]) - (car_ver[0] - ob_rotate_vertex[0][0]) * (ob_rotate_vertex[1][1] - ob_rotate_vertex[0][1])},
                                {(ob_rotate_vertex[2][0] - ob_rotate_vertex[1][0]) * (car_ver[1] - ob_rotate_vertex[1][1]) - (car_ver[0] - ob_rotate_vertex[1][0]) * (ob_rotate_vertex[2][1] - ob_rotate_vertex[1][1])},
                                {(ob_rotate_vertex[3][0] - ob_rotate_vertex[2][0]) * (car_ver[1] - ob_rotate_vertex[2][1]) - (car_ver[0] - ob_rotate_vertex[2][0]) * (ob_rotate_vertex[3][1] - ob_rotate_vertex[2][1])},
                                {(ob_rotate_vertex[0][0] - ob_rotate_vertex[3][0]) * (car_ver[1] - ob_rotate_vertex[3][1]) - (car_ver[0] - ob_rotate_vertex[3][0]) * (ob_rotate_vertex[0][1] - ob_rotate_vertex[3][1])}};
                if(cross_product[0] <= 0.0 && cross_product[1] <= 0.0 && cross_product[2] <= 0.0 && cross_product[3] <= 0.0){
                    return false;
                }
            }
            for(auto& ob_ver : ob_rotate_vertex){
                cross_product = {{(car_rotate_vertex[1][0] - car_rotate_vertex[0][0]) * (ob_ver[1] - car_rotate_vertex[0][1]) - (ob_ver[0] - car_rotate_vertex[0][0]) * (car_rotate_vertex[1][1] - car_rotate_vertex[0][1])},
                                {(car_rotate_vertex[2][0] - car_rotate_vertex[1][0]) * (ob_ver[1] - car_rotate_vertex[1][1]) - (ob_ver[0] - car_rotate_vertex[1][0]) * (car_rotate_vertex[2][1] - car_rotate_vertex[1][1])},
                                {(car_rotate_vertex[3][0] - car_rotate_vertex[2][0]) * (ob_ver[1] - car_rotate_vertex[2][1]) - (ob_ver[0] - car_rotate_vertex[2][0]) * (car_rotate_vertex[3][1] - car_rotate_vertex[2][1])},
                                {(car_rotate_vertex[0][0] - car_rotate_vertex[3][0]) * (ob_ver[1] - car_rotate_vertex[3][1]) - (ob_ver[0] - car_rotate_vertex[3][0]) * (car_rotate_vertex[0][1] - car_rotate_vertex[3][1])}};
                if(cross_product[0] <= 0.0 && cross_product[1] <= 0.0 && cross_product[2] <= 0.0 && cross_product[3] <= 0.0){
                    return false;
                }
            }
        }
    }
    //矩形で衝突判定(低速走行車)
    for(i = 0; i < low_speed_car_states.x.size(); i++){
        low_speed_car_state.x = low_speed_car_states.x[i];
        low_speed_car_state.y = low_speed_car_states.y[i];
        low_speed_car_state.yaw = low_speed_car_states.yaw[i];
        low_speed_car_state.v = low_speed_car_states.v[i];
        low_speed_car_state.rear_x = low_speed_car_states.x[i] - ((WB / 2.f) * cos(low_speed_car_states.yaw[i]));
        low_speed_car_state.rear_y = low_speed_car_states.y[i] - ((WB / 2.f) * sin(low_speed_car_states.yaw[i]));
        for(j = 0; j < fp.x.size(); j++){
            car_vertex = {{fp.x[j] - FULL_LENGTH / 2.f, fp.y[j] - FULL_WIDTH / 2.f},
                        {fp.x[j] - FULL_LENGTH / 2.f, fp.y[j] + FULL_WIDTH / 2.f},
                        {fp.x[j] + FULL_LENGTH / 2.f, fp.y[j] + FULL_WIDTH / 2.f}, 
                        {fp.x[j] + FULL_LENGTH / 2.f, fp.y[j] - FULL_WIDTH / 2.f}};
            car_rotate_vertex = {{(car_vertex[0][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[0][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[0][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[0][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[1][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[1][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[1][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[1][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[2][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[2][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[2][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[2][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]},
                                {(car_vertex[3][0] - fp.x[j]) * cos(fp.yaw[j]) - (car_vertex[3][1] - fp.y[j]) * sin(fp.yaw[j]) + fp.x[j], (car_vertex[3][0] - fp.x[j]) * sin(fp.yaw[j]) + (car_vertex[3][1] - fp.y[j]) * cos(fp.yaw[j]) + fp.y[j]}};
            if(j != 0){
                for(auto k = 0; k < 2; k++){
                    di_target_ind = pure_pursuit_steer_control(low_speed_car_state, target_course, low_speed_car_target_ind);
                    di = di_target_ind.first;
                    low_speed_car_target_ind = di_target_ind.second;
                    low_speed_car_state.update(0.0, di);
                }
            }
            low_speed_car_vertex = {{low_speed_car_state.x - FULL_LENGTH / 2.f, low_speed_car_state.y - FULL_WIDTH / 2.f},
                                    {low_speed_car_state.x - FULL_LENGTH / 2.f, low_speed_car_state.y + FULL_WIDTH / 2.f},
                                    {low_speed_car_state.x + FULL_LENGTH / 2.f, low_speed_car_state.y + FULL_WIDTH / 2.f}, 
                                    {low_speed_car_state.x + FULL_LENGTH / 2.f, low_speed_car_state.y - FULL_WIDTH / 2.f}};
            low_speed_car_rotate_vertex = {{(low_speed_car_vertex[0][0] - low_speed_car_state.x) * cos(low_speed_car_state.yaw) - (low_speed_car_vertex[0][1] - low_speed_car_state.y) * sin(low_speed_car_state.yaw) + low_speed_car_state.x, (low_speed_car_vertex[0][0] - low_speed_car_state.x) * sin(low_speed_car_state.yaw) + (low_speed_car_vertex[0][1] - low_speed_car_state.y) * cos(low_speed_car_state.yaw) + low_speed_car_state.y},
                                            {(low_speed_car_vertex[1][0] - low_speed_car_state.x) * cos(low_speed_car_state.yaw) - (low_speed_car_vertex[1][1] - low_speed_car_state.y) * sin(low_speed_car_state.yaw) + low_speed_car_state.x, (low_speed_car_vertex[1][0] - low_speed_car_state.x) * sin(low_speed_car_state.yaw) + (low_speed_car_vertex[1][1] - low_speed_car_state.y) * cos(low_speed_car_state.yaw) + low_speed_car_state.y},
                                            {(low_speed_car_vertex[2][0] - low_speed_car_state.x) * cos(low_speed_car_state.yaw) - (low_speed_car_vertex[2][1] - low_speed_car_state.y) * sin(low_speed_car_state.yaw) + low_speed_car_state.x, (low_speed_car_vertex[2][0] - low_speed_car_state.x) * sin(low_speed_car_state.yaw) + (low_speed_car_vertex[2][1] - low_speed_car_state.y) * cos(low_speed_car_state.yaw) + low_speed_car_state.y},
                                            {(low_speed_car_vertex[3][0] - low_speed_car_state.x) * cos(low_speed_car_state.yaw) - (low_speed_car_vertex[3][1] - low_speed_car_state.y) * sin(low_speed_car_state.yaw) + low_speed_car_state.x, (low_speed_car_vertex[3][0] - low_speed_car_state.x) * sin(low_speed_car_state.yaw) + (low_speed_car_vertex[3][1] - low_speed_car_state.y) * cos(low_speed_car_state.yaw) + low_speed_car_state.y}};
            for(auto& car_ver : car_rotate_vertex){
                cross_product = {{(low_speed_car_rotate_vertex[1][0] - low_speed_car_rotate_vertex[0][0]) * (car_ver[1] - low_speed_car_rotate_vertex[0][1]) - (car_ver[0] - low_speed_car_rotate_vertex[0][0]) * (low_speed_car_rotate_vertex[1][1] - low_speed_car_rotate_vertex[0][1])},
                                {(low_speed_car_rotate_vertex[2][0] - low_speed_car_rotate_vertex[1][0]) * (car_ver[1] - low_speed_car_rotate_vertex[1][1]) - (car_ver[0] - low_speed_car_rotate_vertex[1][0]) * (low_speed_car_rotate_vertex[2][1] - low_speed_car_rotate_vertex[1][1])},
                                {(low_speed_car_rotate_vertex[3][0] - low_speed_car_rotate_vertex[2][0]) * (car_ver[1] - low_speed_car_rotate_vertex[2][1]) - (car_ver[0] - low_speed_car_rotate_vertex[2][0]) * (low_speed_car_rotate_vertex[3][1] - low_speed_car_rotate_vertex[2][1])},
                                {(low_speed_car_rotate_vertex[0][0] - low_speed_car_rotate_vertex[3][0]) * (car_ver[1] - low_speed_car_rotate_vertex[3][1]) - (car_ver[0] - low_speed_car_rotate_vertex[3][0]) * (low_speed_car_rotate_vertex[0][1] - low_speed_car_rotate_vertex[3][1])}};
                if(cross_product[0] <= 0.0 && cross_product[1] <= 0.0 && cross_product[2] <= 0.0 && cross_product[3] <= 0.0){
                    return false;
                }
            }
            for(auto& low_speed_car_ver : low_speed_car_rotate_vertex){
                cross_product = {{(car_rotate_vertex[1][0] - car_rotate_vertex[0][0]) * (low_speed_car_ver[1] - car_rotate_vertex[0][1]) - (low_speed_car_ver[0] - car_rotate_vertex[0][0]) * (car_rotate_vertex[1][1] - car_rotate_vertex[0][1])},
                                {(car_rotate_vertex[2][0] - car_rotate_vertex[1][0]) * (low_speed_car_ver[1] - car_rotate_vertex[1][1]) - (low_speed_car_ver[0] - car_rotate_vertex[1][0]) * (car_rotate_vertex[2][1] - car_rotate_vertex[1][1])},
                                {(car_rotate_vertex[3][0] - car_rotate_vertex[2][0]) * (low_speed_car_ver[1] - car_rotate_vertex[2][1]) - (low_speed_car_ver[0] - car_rotate_vertex[2][0]) * (car_rotate_vertex[3][1] - car_rotate_vertex[2][1])},
                                {(car_rotate_vertex[0][0] - car_rotate_vertex[3][0]) * (low_speed_car_ver[1] - car_rotate_vertex[3][1]) - (low_speed_car_ver[0] - car_rotate_vertex[3][0]) * (car_rotate_vertex[0][1] - car_rotate_vertex[3][1])}};
                if(cross_product[0] <= 0.0 && cross_product[1] <= 0.0 && cross_product[2] <= 0.0 && cross_product[3] <= 0.0){
                    return false;
                }
            }
        }
    }

    return true;
}

std::vector<FrenetPath> check_paths(std::vector<FrenetPath> fplist, std::vector<std::vector<float>> ob, States low_speed_car_states, TargetCourse target_course, int low_speed_car_target_ind){
    int i;
    int count1 = 0, count2 = 0, count3 = 0, count4 = 0;
    bool ok;
    std::vector<int> ok_ind;

    for(i = 0; i < fplist.size(); i++){
        ok = true;
        for(auto& a : fplist[i].s_dd){
            if(abs(a) > MAX_ACCEL){
                ok = false;
                count1 += 1;
                break;
            }
        }
        if(!ok) continue;
        if(!(check_collision(fplist[i], ob, low_speed_car_states, target_course, low_speed_car_target_ind))){
            count2 += 1;
            continue;
        }
        for(auto& v : fplist[i].s_d){
            if(v > MAX_SPEED){
                ok = false;
                count3 += 1;
                break;
            }
        }
        if(!ok) continue;
        for(auto& c : fplist[i].c){
            if(abs(c) > MAX_CURVATURE){
                ok = false;
                count4 += 1;
                break;
            }
        }
        if(!ok) continue;

        ok_ind.push_back(i);
    }

    // std::cout << count1 << " " << count2 << " " << count3 << " " << count4 << std::endl;

    std::vector<FrenetPath>::iterator fp_it = fplist.begin();
    std::vector<int>::iterator ok_it = ok_ind.begin();
    int fplist_size = fplist.size();
    for(i = 0; i < fplist_size; i++){
        if(i == *ok_it && ok_it < ok_ind.end()){
            ++fp_it;
            ++ok_it;
        }
        else{
            fp_it = fplist.erase(fp_it);
        }
    }

    return fplist;
}

FrenetPath frenet_optimal_planning(TargetCourse target_course, float s0, float c_speed, float c_d, float c_d_d, float c_d_dd, std::vector<std::vector<float>> ob, States low_speed_car_states, int low_speed_car_target_ind){
    std::vector<FrenetPath> fplist;
    float min_cost;
    FrenetPath best_path;

    fplist = calc_frenet_paths(c_speed, c_d, c_d_d, c_d_dd, s0);
    fplist = calc_global_paths(fplist, target_course.csp);
    fplist = check_paths(fplist, ob, low_speed_car_states, target_course, low_speed_car_target_ind);

    // find minimum cost path
    min_cost = std::numeric_limits<float>::infinity();
    for(auto& fp : fplist){
        if(min_cost >= fp.cf){
            min_cost = fp.cf;
            best_path = fp;
        }
    }

    return best_path;
}

TargetCourse::TargetCourse(VectorXd X, VectorXd Y) : csp(X, Y){
    float i_s;
    Vector2d vec;

    for(i_s = 0; i_s < csp.s(csp.s.size() - 1); i_s += 0.1){
        vec = csp.calc_position(i_s);
        rx.push_back(vec(0));
        ry.push_back(vec(1));
        ryaw.push_back(csp.calc_yaw(i_s));
        rk.push_back(csp.calc_curvature(i_s));
    }
    length = 0.0;
    for(auto i = 0; i < rx.size() - 1; i++){
        length += hypot(rx[i] - rx[i + 1], ry[i] - ry[i + 1]);
    }
}

std::pair<int, float> TargetCourse::search_target_index(State state, int old_nearest_point_index){
    std::vector<float> dx, dy, d;
    float distance_this_index, distance_next_index, Lf;
    int i, ind;

    if(old_nearest_point_index == -1){
        for(auto& icx : rx){
            dx.push_back(state.rear_x - icx);
        }
        for(auto& icy : ry){
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
        distance_this_index = state.calc_distance(rx[ind], ry[ind]);
        while(true){
            distance_next_index = state.calc_distance(rx[ind + 1], ry[ind + 1]);
            if(distance_this_index < distance_next_index || std::isnan(distance_next_index)) break;
            if((ind + 1) < rx.size()) ind++;
            distance_this_index = distance_next_index;
        }
        old_nearest_point_index = ind;
    }

    Lf = k * state.v + Lfc;  // update look ahead distance

    while(Lf > state.calc_distance(rx[ind], ry[ind])){
        if((ind + 1) >= rx.size()) break;
        ind++;
    }

    return std::make_pair(ind, Lf);
}

State::State() : x(0.0) ,y(0.0), yaw(0.0), v(0.0), rear_x(x - ((WB /2.f) * cos(yaw))), rear_y(y - ((WB /2.f) * sin(yaw))){}
State::State(float x_, float y_, float yaw_, float v_) : x(x_), y(y_), yaw(yaw_), v(v_), rear_x(x - ((WB /2.f) * cos(yaw))), rear_y(y - ((WB /2.f) * sin(yaw))){}
void State::update(float a, float delta){
    x += v * cos(yaw) * dt;
    y += v * sin(yaw) * dt;
    yaw += v / WB * tan(delta) * dt;
    v += a * dt;
    rear_x = x - ((WB /2.f) * cos(yaw));
    rear_y = y - ((WB /2.f) * sin(yaw));
}
float State::calc_distance(float point_x, float point_y){
    float dx, dy;
    dx = rear_x - point_x;
    dy = rear_y - point_y;
    return hypot(dx, dy);
}

void States::append(State state){
    x.push_back(state.x);
    y.push_back(state.y);
    yaw.push_back(state.yaw);
    v.push_back(state.v);
}
void States::clear(){
    x.clear();
    y.clear();
    yaw.clear();
    v.clear();
}

std::pair<int, float> FrenetPath::search_target_index(State state, int old_nearest_point_index){
    std::vector<float> dx, dy, d;
    float distance_this_index, distance_next_index, Lf;
    int i, ind;

    if(old_nearest_point_index == -1){
        for(auto& icx : x){
            dx.push_back(state.rear_x - icx);
        }
        for(auto& icy : y){
            dy.push_back(state.rear_y - icy);
        }
        for(i = 0; i < dx.size(); i++){
            d.push_back(hypot(dx[i], dy[i]));
        }
        ind = std::min_element(d.begin(),d.end()) - d.begin();
    }
    else{
        ind = old_nearest_point_index;
        distance_this_index = state.calc_distance(x[ind], y[ind]);
        while(true){
            distance_next_index = state.calc_distance(x[ind + 1], y[ind + 1]);
            if(distance_this_index < distance_next_index || std::isnan(distance_next_index)) break;
            if((ind + 1) < x.size()) ind++;
            if((ind + 1) >= x.size()) break;
            distance_this_index = distance_next_index;
        }
    }

    Lf = k * state.v + Lfc;  // update look ahead distance

    while(Lf > state.calc_distance(x[ind], y[ind])){
        if((ind + 1) >= x.size()) break;
        ind++;
    }

    return std::make_pair(ind, Lf);
}


int main(void){
    std::cout << __FILE__ << " start!!" << std::endl;

    int i, j;
    float c_speed, c_d, c_d_d, c_d_dd, s0, area;
    std::vector<float> ob_x, ob_y, low_speed_car_x, low_speed_car_y;
    std::vector<std::vector<float>> ob;
    FrenetPath path;
    // Eigen::VectorXd wx(2), wy(2);
    // wx << 0.0, 50.0;
    // wy << 0.0, 0.0;
    // Eigen::VectorXd wx(5), wy(5);
    // wx << 0.0, 1.0, 2.05, 3.50, 7.05;
    // wy << 0.0, -0.6, 0.5, 0.65, 0.0;
    Eigen::VectorXd wx(9), wy(9);
    wx << 0, 7.696059, 10, 7.380889, 0, -5.843047, -10, -7.071067, 0;
    wy << -10, -6.385191, 0, 6.747034, 10, 8.115343, 0, -7.071067, -10;

    TargetCourse target_course(wx, wy);

    // ob = {{20.0, 10.0}, {30.0, 6.0}, {30.0, 8.0}, {35.0, 8.0}, {50.0, 3.0}};
    // ob = {{30.0, -0.5}, {35.0, -0.5}};
    ob = {{-5.843047, 8.115343}, {-10, 0}};
    // for(float obx = 2; obx < 6; obx += 1){
    //     for(float oby = 0.3; oby < 0.4; oby += 0.1){
    //         ob.push_back({obx, oby});
    //     }
    // }

    for(auto& x : ob){
        ob_x.push_back(x[0]);
        ob_y.push_back(x[1]);
    }

    c_speed = START_SPEED;  // current speed [m/s]
    c_d = 0.0;  // current lateral position [m]
    c_d_d = 0.0;  // current lateral speed [m/s]
    c_d_dd = 0.0;  // current lateral acceleration [m/s]
    s0 = 0.0;  // current course position

    area = 2.0;  // animation area length [m]

    //車の動作
    int car_target_ind = -1, low_speed_car_target_ind, fp_current_ind, target_course_nearest_ind = 0, count = 0, flag = 0;
    float ai, di, targetpath_distance = 0.0, prev_distance, next_distance, cross_product, course_pos = 0.0, min_distance, min_distance1, min_distance2, a, b, c, x1, y1, x2, y2, a1, a2, b1, b2, p1, p2, x, y;
    bool goal = false;
    std::vector<float> ob_distance0, ob_distance1, target_x, target_y, fp_distance, aa, bb;
    std::vector<std::vector<float>> car_vertex, car_rotate_vertex, low_speed_car_vertex, low_speed_car_rotate_vertex;
    std::pair<float, int> di_target_ind;
    std::vector<double> cc,dd;
    Eigen::Vector2d vec;
    std::map<std::string,std::string> keywords, keywords2, keywords3;
    keywords["linewidth"] = "5";
    keywords2["marker"] = "v";
    keywords2["color"] = "b";
    keywords2["markersize"] = "3";
    keywords3["marker"] = "v";
    keywords3["markersize"] = "10";

    State low_speed_car_state1(target_course.rx[OB_POS], target_course.ry[OB_POS], target_course.ryaw[OB_POS], OB_SPEED);
    States car_states, low_speed_car_states;
    low_speed_car_x.push_back(low_speed_car_state1.x);
    low_speed_car_y.push_back(low_speed_car_state1.y);
    low_speed_car_states.append(low_speed_car_state1);

    low_speed_car_target_ind = target_course.search_target_index(low_speed_car_state1, -1).second;
    
    path = frenet_optimal_planning(target_course, s0, c_speed, c_d, c_d_d, c_d_dd, ob, low_speed_car_states, low_speed_car_target_ind);
    
    State car_state(path.x[0], path.y[0], path.yaw[0], c_speed);
    car_states.append(car_state);

    for(i = 0; i < SIM_LOOP; i++){
        std::cout << "path update " << i << std::endl;

        // std::cout << car_state.x << " " << car_state.y << std::endl;
        path = frenet_optimal_planning(target_course, s0, c_speed, c_d, c_d_d, c_d_dd, ob, low_speed_car_states, low_speed_car_target_ind);    
        // std::cout << i << " path.x = " << path.x[0] << " path.y = " << path.y[0] << " s0 = " << s0 << " c_d = " << c_d << std::endl;

        //デバッグ用
        if(i != 0){
            vec = target_course.csp.calc_position(course_pos);
            cc = {vec(0)};
            dd = {vec(1)};
            aa = {car_states.x.back()};
            bb = {car_states.y.back()};
        }

        car_target_ind = path.search_target_index(car_state, car_target_ind).second;

        //UPDATE_PATH_INDまで動いたらpath更新
        // while(UPDATE_PATH_IND > fp_current_ind){
        //UPDATE_STEP(偶数)分動いたらpath更新,pure pursuitはdt=0.1,frenetはDT=0.2
        for(auto j = 0; j < UPDATE_STEP; j++){
            
            ai = proportional_control(path.s_d[car_target_ind], car_state.v);

            //car
            di_target_ind = pure_pursuit_steer_control(car_state, path, car_target_ind);
            di = di_target_ind.first;
            car_target_ind = di_target_ind.second;
            car_state.update(ai, di);  // Control vehicle
            target_x = {path.x[car_target_ind]};
            target_y = {path.y[car_target_ind]};

            //low speed car
            di_target_ind = pure_pursuit_steer_control(low_speed_car_state1, target_course, low_speed_car_target_ind);
            di = di_target_ind.first;
            low_speed_car_target_ind = di_target_ind.second;
            low_speed_car_state1.update(0.0, di);
            low_speed_car_x = {low_speed_car_state1.x};
            low_speed_car_y = {low_speed_car_state1.y};

            fp_distance.clear();
            for(auto ind = 0; ind < path.x.size(); ind++){
                fp_distance.push_back(hypot(car_state.x - path.x[ind], car_state.y - path.y[ind]));
            }
            fp_current_ind = std::min_element(fp_distance.begin(),fp_distance.end()) - fp_distance.begin();

            car_states.append(car_state);
            low_speed_car_states.clear();
            low_speed_car_states.append(low_speed_car_state1);

            // 障害物との距離での評価
            // ob_distance0.push_back(hypot(ob_x[0] - car_state.x, ob_y[0] - car_state.y));
            // ob_distance1.push_back(hypot(ob_x[1] - car_state.x, ob_y[1] - car_state.y));

            // 目標パスとの距離での評価
            // if(car_state.y < 0){
            //     targetpath_distance += (car_state.y * -1);
            // }
            // else{
            //     targetpath_distance += car_state.y;
            // }
            // count++;

            if(show_animation){
                plt::cla();
                plt::plot(target_course.rx, target_course.ry);
                plt::plot(ob_x, ob_y, "xk");
                plt::plot(low_speed_car_x, low_speed_car_y, "s");
                plt::plot(aa, bb, keywords3);
                plt::plot(cc, dd, keywords3);
                plt::plot(path.x, path.y, "-or");
                plt::plot(target_x, target_y, keywords2);
                plt::plot(car_states.x, car_states.y, keywords);
                plt::xlim(car_state.x - area, car_state.x + area);
                plt::ylim(car_state.y - area, car_state.y + area);
                plt::xlabel("x[m]");
                plt::ylabel("y[m]");
                plt::grid(true);
                plt::pause(0.0001);
            }

            if(hypot(car_state.x - target_course.rx.back(), car_state.y - target_course.ry.back()) <= ROBOT_RADIUS && course_pos >= target_course.length * 0.9){
                std::cout << "Goal" << std::endl;
                goal = true;
                // std::ofstream ofs("csv/experiment4.csv", std::ios::app);
                //         for(auto& obx : car_states.x){
                //             ofs << obx << ",";
                //         }
                //         ofs << std::endl;
                //         for(auto& oby : car_states.y){
                //             ofs << oby << ",";
                //         }
                //         ofs << std::endl;
                // std::ofstream ofs("../PathPlanning/ob_distance.csv", std::ios::app);
                // ofs << *std::min_element(ob_distance0.begin(), ob_distance0.end()) << "," << *std::min_element(ob_distance1.begin(), ob_distance1.end()) << std::endl;
                // std::ofstream o("../PathPlanning/targetpath_distance.csv", std::ios::app);
                // o << START_SPEED * 3.6 << ", " << TARGET_SPEED * 3.6 << "," << targetpath_distance / count << std::endl;
                break;
            }
        }

        if(goal) break;
        
        s0 = path.s[UPDATE_STEP / 2];
        c_d = path.d[UPDATE_STEP / 2];
        c_d_d = path.d_d[UPDATE_STEP / 2];
        c_d_dd = path.d_dd[UPDATE_STEP / 2];
        c_speed = path.s_d[UPDATE_STEP / 2];
        Vector2d vect = target_course.csp.calc_position(s0);
        // std::cout << "x = " << vect(0) << " y = " << vect(1) << " s0 = " << s0 << std::endl;

        for(auto j = target_course_nearest_ind; j < target_course.rx.size() - 1; j++){
            prev_distance = hypot(car_state.x - target_course.rx[j], car_state.y - target_course.ry[j]);
            next_distance = hypot(car_state.x - target_course.rx[j + 1], car_state.y - target_course.ry[j + 1]);
            if(prev_distance < next_distance){
                target_course_nearest_ind = j;
                a1 = target_course.rx[target_course_nearest_ind];
                a2 = target_course.ry[target_course_nearest_ind];
                b1 = target_course.rx[target_course_nearest_ind + 1];
                b2 = target_course.ry[target_course_nearest_ind + 1];
                p1 = car_state.x;
                p2 = car_state.y;
                a = a2 - b2;
                b = b1 - a1;
                c = a1 * (b2 - a2) - a2 * (b1 - a1);
                min_distance = (abs(a * p1 + b * p2 + c) / sqrt(pow(a, 2) + pow(b, 2)));
                cross_product = (target_course.rx[target_course_nearest_ind + 1] - target_course.rx[target_course_nearest_ind]) * (car_state.y - target_course.ry[target_course_nearest_ind]) - (car_state.x - target_course.rx[target_course_nearest_ind]) * (target_course.ry[target_course_nearest_ind + 1] - target_course.ry[target_course_nearest_ind]);
                if(cross_product < 0){
                    min_distance *= -1;
                }
                break;
            }
        }
        course_pos = 0.0;
        for(auto j = 0; j < target_course_nearest_ind; j++){
            course_pos += hypot(target_course.rx[j] - target_course.rx[j + 1], target_course.ry[j] - target_course.ry[j + 1]);
        }
        course_pos = course_pos / target_course.length * float(target_course.csp.s(target_course.csp.s.size() - 1));
        
        // std::cout << target_course.rx.size() << " / " << target_course_nearest_ind << " x = " << car_state.x << " y = " << car_state.y << " s0 = " << course_pos << std::endl;
        // if(hypot(car_state.x - path.x[car_target_ind], car_state.y - path.y[car_target_ind]) < 0.3 || hypot(car_state.x - path.x[car_target_ind], car_state.y - path.y[car_target_ind]) > 0.5){
            s0 = course_pos;
            c_d = min_distance;
            c_speed = car_state.v;
            // c_d_d = 
            // c_d_dd = 
            // std::cout << s0 << " " << c_d << " " << c_speed << std::endl;
        // }
        // std::cout << path.s[UPDATE_STEP/2] << s0 << std::endl;
    }

    std::cout << "Finish" << std::endl;
    if(show_animation){  // pragma: no cover
        plt::grid(true);
        plt::pause(0.0001);
        plt::show();
    }

    return 0;
}