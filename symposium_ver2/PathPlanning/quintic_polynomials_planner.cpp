#include <iostream>
#include <eigen3/Eigen/Dense>
#include "quintic_polynomials_planner.h"
using namespace std;
using namespace Eigen;

// float MAX_T = 100.0;
// float MIN_T = 5.0;

QuinticPolynomial::QuinticPolynomial(){}
QuinticPolynomial::QuinticPolynomial(float xs, float vxs, float axs, float xe, float vxe, float axe, float time){
    a0 = xs;
    a1 = vxs;
    a2 = axs / 2.0;
    
    Matrix3d A;
    Vector3d b;
    A << pow(time,3), pow(time,4), pow(time,5), 3 * pow(time,2), 4 * pow(time,3), 5 * pow(time,4), 6 * time, 12 * pow(time,2), 20 * pow(time,3);
    b << xe - a0 - a1 * time - a2 * pow(time,2), vxe - a1 - 2 * a2 * time, axe - 2 * a2;
    Vector3d x = A.colPivHouseholderQr().solve(b);
    
    a3 = x[0];
    a4 = x[1];
    a5 = x[2];
}

float QuinticPolynomial::calc_point(float t){
    float xt;
    xt = a0 + a1 * t + a2 * pow(t,2) + a3 * pow(t,3) + a4 * pow(t,4) + a5 * pow(t,5);
    
    return xt;
}

float QuinticPolynomial::calc_first_derivative(float t){
    float xt;
    xt = a1 + 2 * a2 * t + 3 * a3 * pow(t,2) + 4 * a4 * pow(t,3) + 5 * a5 * pow(t,4);

    return xt;
}

float QuinticPolynomial::calc_second_derivative(float t){
    float xt;
    xt = 2 * a2 + 6 * a3 * t + 12 * a4 * pow(t,2) + 20 * a5 * pow(t,3);

    return xt;
}

float QuinticPolynomial::calc_third_derivative(float t){
    float xt;
    xt = 6 * a3 + 24 * a4 * t + 60 * a5 * pow(t,2);

    return xt;
}


// void QuinticPolynomial::quintic_polynomials_planner(float sx, float sy, float syaw, float sv, float sa, float gx, float gy, float gyaw, float gv, float ga, float max_accel, float max_jerk, float dt){
//     float vxs = sv * cos(syaw);
//     float vys = sv * sin(syaw);
//     float vxg = gv * cos(gyaw);
//     float vyg = gv * sin(gyaw);
//     float axs = sa * cos(syaw);
//     float ays = sa * sin(syaw);
//     float axg = ga * cos(gyaw);
//     float ayg = ga * sin(gyaw);

//     int i, k, p = 0;
//     int num = int((MAX_T + dt) / dt);
//     float T, t, vx, vy, v, yaw, ax, ay, a, jx, jy, j;
//     float time_[num];
//     float rx_[num];
//     float ry_[num];
//     float ryaw_[num];
//     float rv_[num];
//     float ra_[num];
//     float rj_[num];

//     for(T = MIN_T; T < MAX_T; T + MIN_T){
//         QuinticPolynomial xqp(sx, vxs, axs, gx, vxg, axg, T);
//         QuinticPolynomial yqp(sx, vxs, axs, gx, vxg, axg, T);

//         for(i = 0; i < num; i++){
//             time_[i] = 0;
//             rx_[i] = 0;
//             ry_[i] = 0;
//             ryaw_[i] = 0;
//             rv_[i] = 0;
//             ra_[i] = 0;
//             rj_[i] = 0;
//         }

//         for(t = 0.0; t < T + dt; t + dt){
//             time_[p] = t;
//             rx_[p] = xqp.calc_point(t);
//             ry_[p] = yqp.calc_point(t);

//             vx = xqp.calc_first_derivative(t);
//             vy = yqp.calc_first_derivative(t);
//             v = hypot(vx, vy);
//             yaw = atan2(vy, vx);
//             rv_[p] = v;
//             ryaw_[p] = yaw;

//             ax = xqp.calc_second_derivative(t);
//             ay = yqp.calc_second_derivative(t);
//             a = hypot(ax, ay);
//             // len(rv) >= 2の書き換え
//             if(rv_[0] != 0 && rv_[1] != 0 && rv_[num-1] - rv_[num-2] < 0.0){
//                 a *= -1;
//             }
//             ra_[p] = a;

//             jx = xqp.calc_third_derivative(t);
//             jy = yqp.calc_third_derivative(t);
//             j = hypot(jx, jy);
//             // len(ra) >= 2の書き換え
//             if(ra_[0] != 0 && ra_[1] != 0 && ra_[num-1] - ra_[num-2] < 0.0){
//                 j *= -1;
//             }
//             rj_[p] = j;

//             p++;
//         }

//         float max_ra = 0;
//         float max_rj = 0;
//         for(i = 0; i < num; i++){
//             if(abs(ra_[i]) > max_ra){
//                 max_ra = abs(ra_[i]);
//             }
//             if(abs(rj_[i]) > max_rj){
//                 max_rj = abs(rj_[i]);
//             }
//         }
//         if(max_ra <= max_accel && max_rj <= max_jerk){
//             cout << "find path!!" << endl;
//             break;
//         }
//     }

//     time = time_;
//     rx = rx_;
//     ry = ry_;
//     ryaw = ryaw_;
//     rv = rv_;
//     ra = ra_;
//     rj = rj_;
// }

// int main(void){
//     float c_d;  // current lateral position [m]
//     float c_d_d;  // current lateral speed [m/s]
//     float c_d_dd;  // current lateral acceleration [m/s]
//     float di;
//     float Ti;

//     int i = 0;
//     for(i = 0; i < 30; i++){
//         cin >> c_d >> c_d_d >> c_d_dd >> di >> Ti;

//         QuinticPolynomial quintic_polynomial(c_d, c_d_d, c_d_dd, di, 0.0, 0.0, Ti);
//         cout << quintic_polynomial.a0 << " " << quintic_polynomial.a1 << " " << quintic_polynomial.a2 << " " << quintic_polynomial.a3 << " " << quintic_polynomial.a4 << " " << quintic_polynomial.a5 << endl;
//     }

//     return 0;
// }