#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include <Python.h>
#include "matplotlibcpp.h"
#include <string>
#include <algorithm>

namespace plt = matplotlibcpp;

// Generate 2D parabola
double fx(double x){
    return pow(x, 2);
}

// Generate 3D parabola
double fx3D(double x, double y){
    return pow(x, 2) + pow(y, 2);
}

// Build a grid map of the original 3D parabola equation
std::map<std::string, std::vector<std::vector<double>>> arange3D(double a, double b){
    std::map<std::string, std::vector<std::vector<double>>> z;
    int n = 50;
    double dx = (b - a)/((double) n - 1);
    std::vector<double> tx, ty, tz;
    for(int i = 0; i < n; ++i){
        tx.clear();
        ty.clear();
        tz.clear();
        for(int j = 0; j < n; ++j){
            tx.push_back(a + i*dx);
            ty.push_back(a + j*dx);
            tz.push_back(fx3D(a + i*dx, a + j*dx));
        }
        z["x"].push_back(tx);
        z["y"].push_back(ty);
        z["z"].push_back(tz);
    }
    return z;
}

// Build a line between the points a and b
std::vector<double> arange(double a, double b){
    int n = 40;
    double dx = (b - a)/((double) n - 1);
    std::vector<double> result;
    for(int i = 0; i < n; ++i){
        result.push_back(a + i*dx);
    }
    return result;
}

// Matrix multiplication function
std::vector<std::vector<double>> mmult(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y){
    std::vector<std::vector<double>> z;
    std::vector<double> temp;
    double total = 0;
    for(int i = 0; i < x.size(); ++i){
        temp.clear();
        for(int j = 0; j < y[0].size(); ++j){
            total = 0;
            for(int k = 0; k < x[0].size(); ++k){
                total += x[i][k]*y[k][j];
            }
            temp.push_back(total);
        }
        z.push_back(temp);
    }
    return z;
}

// Rotation matrix which uses sine and cosine to rotate by an inputted theta variable
void Rotate(std::vector<double> & x, std::vector<double> & y, double theta){
    for(int i = 0; i < x.size(); ++i){
        std::vector<std::vector<double>> temp, rf;
        rf = {{std::cos(theta), -std::sin(theta)},{std::sin(theta), std::cos(theta)}};
        temp = {{x[i]}, {y[i]}};
        temp = mmult(rf, temp);
        x[i] = temp[0][0];
        y[i] = temp[1][0];
    }
}

// 3D Rotation matrix which uses sine and cosine to rotate 3D function with an inputted theta
void Rotate3D(std::map<std::string, std::vector<std::vector<double>>> & H, double theta){
    std::vector<std::vector<double>> rotateX, rotateY, rotateZ, turn;

    // Rotate along x-axis
    rotateX = {
        {1, 0, 0},
        {0, std::cos(theta), -std::sin(theta)},
        {0, std::sin(theta), std::cos(theta)}
    };

    // Rotate along y-axis
    rotateY = {
        {std::cos(theta), 0, -std::sin(theta)},
        {0, 1, 0},
        {std::sin(theta), 0, std::cos(theta)}
    };

    // Rotate along z-axis
    rotateZ = {
        {std::cos(theta), -std::sin(theta), 0},
        {std::sin(theta), std::cos(theta), 0},
        {0, 0, 1}
    };

    for(int i = 0; i < H["x"].size(); ++i){
        for(int j = 0; j < H["x"][0].size(); ++j){
            // Rotate each point
            turn = {{H["x"][i][j]}, {H["y"][i][j]}, {H["z"][i][j]}};
            turn = mmult(rotateX, turn);
            turn = mmult(rotateY, turn);
            turn = mmult(rotateZ, turn);

            // Update each point after rotation
            H["x"][i][j] = turn[0][0];
            H["y"][i][j] = turn[1][0];
            H["z"][i][j] = turn[2][0];
        }
    }

}

// 2D rotation animation
void rotate2D()
{
    // Declare 2D plot
    PyObject * ax = plt::chart2D(111);

    std::vector<double> x, y;

    // Set bounds
    x = arange(-4, 4);
    for(auto & i : x){
        y.push_back(fx(i));
    }

    // Animated rotation plot
    for(int i = 0; i < 20; ++i){
        plt::Clear3DChart(ax);
        Rotate(std::ref(x), std::ref(y), 0.05);
        plt::plot2D(ax, x, y, "red");
        plt::pause(1);
    }

    plt::show();
    
}

// 3D rotation animation
void rotate3D()
{
    // Generate grid map
    std::map<std::string, std::vector<std::vector<double>>> H = arange3D(-4, 4);

    // Initialize 3D plot
    PyObject * ax = plt::chart(111);

    // Animated plot
    for(int i = 0; i < 40; ++i){
        plt::Clear3DChart(ax);

        // Store data from each rotation
        Rotate3D(std::ref(H), 0.05);

        // Plot the rotation
        plt::surface3DMap(ax, H["x"], H["y"], H["z"], "jet", 1.0);
        plt::pause(1);
    }

    plt::show();
}


int main()
{
    // Animate 3D rotation plot
    rotate3D();

    return 0;
}
