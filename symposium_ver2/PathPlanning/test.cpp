#include <iostream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h>


int main(void){
    float FULL_LENGTH = 4.8;
    float FULL_WIDTH = 1.7;
    float theta = -M_PI / 2;
    std::vector<float> center;
    std::vector<std::vector<float>> vertex, rotate_vertex;
    center = {1,1};
    vertex = {{center[0] - FULL_LENGTH / 2, center[1] - FULL_WIDTH / 2}, 
            {center[0] - FULL_LENGTH / 2, center[1] + FULL_WIDTH / 2}, 
            {center[0] + FULL_LENGTH / 2, center[1] + FULL_WIDTH / 2}, 
            {center[0] + FULL_LENGTH / 2, center[1] - FULL_WIDTH / 2}};

    rotate_vertex = {{(vertex[0][0] - center[0]) * std::cos(theta) - (vertex[0][1] - center[1]) * std::sin(theta) + center[0], (vertex[0][0] - center[0]) * std::sin(theta) + (vertex[0][1] - center[1]) * std::cos(theta) + center[1]}, 
                    {(vertex[1][0] - center[0]) * std::cos(theta) - (vertex[1][1] - center[1]) * std::sin(theta) + center[0], (vertex[1][0] - center[0]) * std::sin(theta) + (vertex[1][1] - center[1]) * std::cos(theta) + center[1]}, 
                    {(vertex[2][0] - center[0]) * std::cos(theta) - (vertex[2][1] - center[1]) * std::sin(theta) + center[0], (vertex[2][0] - center[0]) * std::sin(theta) + (vertex[2][1] - center[1]) * std::cos(theta) + center[1]}, 
                    {(vertex[3][0] - center[0]) * std::cos(theta) - (vertex[3][1] - center[1]) * std::sin(theta) + center[0], (vertex[3][0] - center[0]) * std::sin(theta) + (vertex[3][1] - center[1]) * std::cos(theta) + center[1]}};

    for(auto& i : rotate_vertex){
        std::cout << i[0] << " " << i[1] << std::endl;
    }

    return 0;
}