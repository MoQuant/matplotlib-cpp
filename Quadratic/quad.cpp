#define USE_MATH_DEFINES
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <Python.h>
#include "matplotlibcpp.h"


namespace plt = matplotlibcpp;

double pi(){
    return M_PI;
}

double FX3D(double x, double y){
    return std::sin(x*y) + std::cos(x*y);
}

std::vector<std::vector<double>> MMULT(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y){
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

std::vector<std::vector<double>> Matrix(double a, double b){
    return {{a}, {b}};
}

std::vector<std::vector<double>> TRANSPOSE(std::vector<std::vector<double>> x){
    std::vector<std::vector<double>> y;
    std::vector<double> ty;
    for(int i = 0; i < x[0].size(); ++i){
        ty.clear();
        for(int j = 0; j < x.size(); ++j){
            ty.push_back(x[j][i]);
        }
        y.push_back(ty);
    }
    return y;
}

std::vector<std::vector<double>> Jacobian(double x0, double y0){
    double first = std::cos(x0*y0)*y0;
    double second = std::sin(x0*y0)*y0;
    double third = std::cos(x0*y0)*x0;
    double fourth = std::sin(x0*y0)*x0;
    return {{first - second}, {third - fourth}};
}

std::vector<std::vector<double>> Hessian(double x0, double y0){
    double first = -std::sin(x0*y0)*pow(y0, 2) - std::cos(x0*y0)*pow(y0, 2);
    double second = -std::sin(x0*y0)*pow(x0, 2) - std::cos(x0*y0)*pow(x0, 2);
    double third = std::cos(x0*y0) - std::sin(x0*y0)*x0*y0 - std::sin(x0*y0) - std::cos(x0*y0)*x0*y0;
    return {{first, third},{third, second}};
}

std::map<std::string, std::vector<std::vector<double>>> QUADRATIC(std::map<std::string, std::vector<std::vector<double>>> EQ, double x0, double y0){
    std::map<std::string, std::vector<std::vector<double>>> QFIN;
    std::vector<std::vector<double>> D2;

    std::vector<double> tx, ty, tz;
    for(int i = 0; i < EQ["x"].size(); ++i){
        tx.clear();
        ty.clear();
        tz.clear();
        for(int j = 0; j < EQ["x"][0].size(); ++j){
            D2 = Matrix(EQ["x"][i][j] - x0, EQ["y"][i][j] - y0);
            double atime = FX3D(x0, y0);
            double jtime = MMULT(TRANSPOSE(Jacobian(x0, y0)), D2)[0][0];
            double htime = MMULT(TRANSPOSE(D2), MMULT(Hessian(x0, y0), D2))[0][0];
            tx.push_back(EQ["x"][i][j]);
            ty.push_back(EQ["y"][i][j]);
            tz.push_back(atime + jtime + htime);
        }
        QFIN["x"].push_back(tx);
        QFIN["y"].push_back(ty);
        QFIN["z"].push_back(tz);
    }

    return QFIN;
}

std::map<std::string, std::vector<std::vector<double>>> GRID(double a, double b){
    std::map<std::string, std::vector<std::vector<double>>> G;
    int n = 50;
    double dx = (b - a)/((double) n - 1);
    std::vector<double> tx, ty, tz;
    for(int i = 0; i < n; ++i){
        double u = a + i*dx;
        tx.clear();
        ty.clear();
        tz.clear();
        for(int j = 0; j < n; ++j){
            double v = a + j*dx;
            tx.push_back(u);
            ty.push_back(v);
            tz.push_back(FX3D(u, v));
        }
        G["x"].push_back(tx);
        G["y"].push_back(ty);
        G["z"].push_back(tz);
    }
    return G;
}

int main()
{
    PyObject * ax = plt::chart(111);

    std::map<std::string, std::vector<std::vector<double>>> EQ = GRID(-pi()/2, pi()/2);
    std::map<std::string, std::vector<std::vector<double>>> QuantFinance = QUADRATIC(EQ, 1.0, 1.0);

    plt::surface3DMap(ax, EQ["x"], EQ["y"], EQ["z"], "jet", 1.0);
    plt::surface3D(ax, QuantFinance["x"], QuantFinance["y"], QuantFinance["z"], "red", 1.0);

    plt::show();

    return 0;
}