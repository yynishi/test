#include <iostream>
#include <eigen3/Eigen/Dense>
#include "cubic_spline_planner.h"
using namespace std;
using namespace Eigen;

Spline::Spline(VectorXd X, VectorXd Y){
    int i;
    float tb;
    MatrixXd A;
    nx = X.rows();
    VectorXd h(nx - 1), tmp_b(nx - 1), tmp_d(nx - 1), B;
    
    x = X;
    y = Y;

    for(i = 0; i < nx - 1; i++){
        h(i) = x(i + 1) - x(i);
    }

    A = this->__calc_A(h);
    B = this->__calc_B(h);
    c = A.colPivHouseholderQr().solve(B);

    for(i = 0; i < nx - 1; i++){
        tmp_d(i) = (c(i + 1) - c(i)) / (3.0 * h(i)); 
        tb = (y(i + 1) - y(i)) / h(i) - h(i) * (c(i + 1) + 2.0 * c(i)) / 3.0;
        tmp_b(i) = tb;
    }
    b = tmp_b;
    d = tmp_d;
}

MatrixXd Spline::__calc_A(VectorXd h){
    int i;
    MatrixXd A(nx, nx);

    A(0,0) = 1.0;
    for(i = 0; i < nx - 1; i++){
        if(i != nx - 2){
            A(i + 1, i + 1) = 2.0 * (h(i) + h(i + 1));
        }
        A(i + 1, i) = h(i);
        A(i, i + 1) = h(i);
    }

    A(0, 1) = 0.0;
    A(nx - 1, nx - 2) = 0.0;
    A(nx - 1, nx - 1) = 1.0;

    return A;
}

VectorXd Spline::__calc_B(VectorXd h){
    int i;
    VectorXd B(nx);

    for(i = 0; i < nx - 2; i++){
        B(i + 1) = 3.0 * (y(i + 2) - y(i + 1)) / h(i + 1) - 3.0 * (y(i + 1) - y(i)) / h(i);
    }

    return B;
}

float Spline::calc(float t){
    if(t < x(0))
        return NAN;
    else if(t > x(nx - 1))
        return NAN;
    
    int i;
    float dx, result;

    i = __search_index(t);
    dx = t - x(i);
    result = y(i) + b(i) * dx + c(i) * pow(dx, 2) + d(i) * pow(dx, 3);

    return result;
}

float Spline::calcd(float t){
    if(t < x(0))
        return NAN;
    else if(t > x(nx - 1))
        return NAN;
    
    int i;
    float dx, result;

    i = __search_index(t);
    dx = t - x(i);
    result = b(i) + 2.0 * c(i) * dx + 3.0 * d(i) * pow(dx, 2);

    return result;
}

float Spline::calcdd(float t){
    if(t < x(0))
        return NAN;
    else if(t > x(nx - 1))
        return NAN;
    
    int i;
    float dx, result;

    i = __search_index(t);
    dx = t - x(i);
    result = 2.0 * c(i) + 6.0 * d(i) * dx;

    return result;
}

int Spline::__search_index(float a){
    int i, index, len = x.rows();

    for(i = 0; i < len; i++){
        if(x(i) <= a){
            index = i;
        }
    }

    return index;
}


Spline2D::Spline2D(VectorXd X, VectorXd Y) : s(__calc_s(X, Y)), sx(s, X), sy(s, Y){}

VectorXd Spline2D::__calc_s(VectorXd X, VectorXd Y){
    int i, nx = X.rows();
    float total = 0, dx, dy, ds;
    VectorXd S(nx);
    
    S(0) = 0;
    for(i = 0; i < nx - 1; i++){
        dx = X(i + 1) - X(i);
        dy = Y(i + 1) - Y(i);
        total += hypot(dx, dy);
        S(i + 1) = total;
    }

    return S;
}

Vector2d Spline2D::calc_position(float S){
    Vector2d vec;
    vec(0) = sx.calc(S);
    vec(1) = sy.calc(S);
    return vec;
}

float Spline2D::calc_curvature(float S){
    float dx, ddx, dy, ddy, k;
    dx = sx.calcd(S);
    ddx = sx.calcdd(S);
    dy = sy.calcd(S);
    ddy = sy.calcdd(S);
    k = (ddy * dx - ddx * dy) / (pow(pow(dx, 2) + pow(dy, 2), 1.5));
    return k;
}

float Spline2D::calc_yaw(float S){
    float dx, dy, yaw;
    dx = sx.calcd(S);
    dy = sy.calcd(S);
    yaw = atan2(dy, dx);
    return yaw;
}

// int main(void){
//     float i;
//     VectorXd X(5), Y(5);
//     float yaw;
//     X << 0.0, 10.0, 20.5, 35.0, 70.5;
//     Y << 0.0, -6.0, 5.0, 6.5, 0.0;
//     Spline2D csp(X,Y);
//     cout << csp.s << endl;
// }