#ifndef _QUINTIC_POLYNOMIALS_PLANNER_H_
#define _QUINTIC_POLYNOMIALS_PLANNER_H_

class QuinticPolynomial{
    public:
        float a0;
        float a1;
        float a2;
        float a3;
        float a4;
        float a5;

        // float *time = NULL;
        // float *rx = NULL;
        // float *ry = NULL;
        // float *ryaw = NULL;
        // float *rv = NULL;
        // float *ra = NULL;
        // float *rj = NULL;
    public:
        QuinticPolynomial();
        QuinticPolynomial(float xs, float vxs, float axs, float xe, float vxe, float axe, float time);
        float calc_point(float t);
        float calc_first_derivative(float t);
        float calc_second_derivative(float t);
        float calc_third_derivative(float t);
        // void quintic_polynomials_planner(float sx, float sy, float syaw, float sv, float sa, float gx, float gy, float gyaw, float gv, float ga, float max_accel, float max_jerk, float dt);
};


#endif // _QUINTIC_POLYNOMIALS_PLANNER_H_