#include <iostream>
#include <vector>
#include <string>
#include <Python.h>
#include "matplotlibcpp.h"
#include <math.h>
#include <map>

namespace plt = matplotlibcpp;

// Parabola Function
double fx(double x, double y)
{
    return -(pow(x, 2) + pow(y, 2));
}

// Gradient for dx
double dfx(double x, double y)
{
    return -2.0*x;
}

// Gradient for dy
double dfy(double x, double y)
{
    return -2.0*y;
}

// Plane Equation
double plane(double x, double y, double x0, double y0)
{
    return dfx(x0, y0)*(x - x0) + dfy(x0, y0)*(y - y0) + fx(x0, y0);
}

// Generates a tangent plane
std::map<std::string, std::vector<std::vector<double>>> Tangent(double c, double x0, double y0)
{
    // Setting grid bounds
    double m0 = x0 - c, m1 = x0 + c;
    double n0 = y0 - c, n1 = y0 + c;
    int o = 60;

    double dM = (m1 - m0)/(o - 1);
    double dN = (n1 - n0)/(o - 1);

    std::map<std::string, std::vector<std::vector<double>>> result;
    std::vector<double> tx, ty, tz;

    double rx, ry;

    // Building the plane at a centered point (where the current gradient is)
    for(int i = 0; i < o; ++i){
        tx.clear();
        ty.clear();
        tz.clear();
        rx = m0 + i*dM;
        for(int j = 0; j < o; ++j){
            ry = n0 + j*dM;
            tx.push_back(rx);
            ty.push_back(ry);
            tz.push_back(plane(rx, ry, x0, y0));
        }
        result["x"].push_back(tx);
        result["y"].push_back(ty);
        result["z"].push_back(tz);
    }
    return result;
}

int main()
{
    // Declare the plot
    PyObject * ax = plt::chart(111);
    plt::Clear3DChart(ax);

    // Define range
    int n = 60;
    double t0 = -4.0, t1 = 4.0;
    double dT = (t1 - t0)/(n - 1);

    std::vector<std::vector<double>> x, y, z;
    std::vector<double> tx, ty, tz;

    double rx, ry;

    // Build 3D graph for the parabola
    for(int i = 0; i < n; ++i){
        tx.clear();
        ty.clear();
        tz.clear();
        rx = t0 + i*dT;
        for(int j = 0; j < n; ++j){
            ry = t0 + j*dT;
            tx.push_back(rx);
            ty.push_back(ry);
            tz.push_back(fx(rx, ry));
        }
        x.push_back(tx);
        y.push_back(ty);
        z.push_back(tz);
    }

    // Build 3D animation for the gradient reaching its maximum and the plane centered around the gradient points
    double Px = -3.5, Py = -3.5;
    double learning = 0.15;
    
    std::map<std::string, std::vector<std::vector<double>>> gd;

    for(int i = 0; i < 100; ++i){
        plt::Clear3DChart(ax);
        gd = Tangent(2, Px, Py);

        // Plots the parabola and gradient plane
        plt::surface3D(ax, x, y, z, "red", 0.9);
        plt::surface3D(ax, gd["x"], gd["y"], gd["z"], "green", 0.9);

        // Updates the gradients in gradient descent
        Px = Px + learning*dfx(Px, Py);
        Py = Py + learning*dfy(Px, Py);

        plt::pause(0.5);
    }


    plt::show();

    return 0;
}
