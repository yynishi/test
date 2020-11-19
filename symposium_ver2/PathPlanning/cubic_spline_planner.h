#ifndef _CUBIC_SPLINE_PLANNER_H_
#define _CUBIC_SPLINE_PLANNER_H_

#include <eigen3/Eigen/Dense>
using namespace Eigen;

class Spline{
    public:
        VectorXd x, y, b, c, d, w;
        int nx;
    public:
        Spline(VectorXd X, VectorXd Y);
        float calc(float t);
        float calcd(float t);
        float calcdd(float t);
        int __search_index(float a);
        MatrixXd __calc_A(VectorXd h);
        VectorXd __calc_B(VectorXd h);
};

class Spline2D{
    public:
        VectorXd s;
        Spline sx, sy;
    public:
        Spline2D(VectorXd X, VectorXd Y);
        VectorXd __calc_s(VectorXd X, VectorXd Y);
        Vector2d calc_position(float S);
        float calc_curvature(float S);
        float calc_yaw(float S);
};

#endif // _CUBIC_SPLINE_PLANNER_H_