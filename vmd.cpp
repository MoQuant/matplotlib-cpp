#include <iostream>
#include <vector>
#include <Python.h>
#include "matplotlibcpp.h"
#include <math.h>
#include <cmath>

namespace plt = matplotlibcpp;

void PRINTX(std::vector<std::vector<double>> x){
    for(auto & row : x){
        for(auto & col : row){
            std::cout << col << "\t";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> MMULT(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y){
    std::vector<std::vector<double>> z;
    std::vector<double> t;
    double total = 0;
    for(int i = 0; i < x.size(); ++i){
        t.clear();
        for(int j = 0; j < y[0].size(); ++j){
            total = 0;
            for(int k = 0; k < x[0].size(); ++k){
                total += x[i][k]*y[k][j];
            }
            t.push_back(total);
        }
        z.push_back(t);
    }
    return z;
}

std::vector<std::vector<double>> INVERSE(std::vector<std::vector<double>> x){
    std::vector<std::vector<double>> I;
    std::vector<double> temp;
    int n = x.size();
    for(int i = 0; i < n; ++i){
        temp.clear();
        for(int j = 0; j < n; ++j){
            if(i == j){
                temp.push_back(1.0);
            } else {
                temp.push_back(0.0);
            }
        }
        I.push_back(temp);
    }

    for(int i = 1; i < n; ++i){
        for(int j = 0; j < i; ++j){
            double A = x[i][j];
            double B = x[j][j];
            for(int k = 0; k < n; ++k){
                x[i][k] -= (A/B)*x[j][k];
                I[i][k] -= (A/B)*I[j][k];
            }
        }
    }

    for(int i = 1; i < n; ++i){
        for(int j = 0; j < i; ++j){
            double A = x[j][i];
            double B = x[i][i];
            for(int k = 0; k < n; ++k){
                x[j][k] -= (A/B)*x[i][k];
                I[j][k] -= (A/B)*I[i][k];
            }
        }
    }

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            I[i][j] /= x[i][i];
        }
    }
    
    return I;
}

int main()
{
    PyObject * ax = plt::chart2D(111);

    std::vector<double> x = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {5, 6, 3, 1, 2, 9, 1, 2, 5, 7};

    std::vector<std::vector<double>> X, Y;
    std::vector<double> temp;

    for(int i = 0; i < x.size(); ++i){
        temp.clear();
        for(int j = 0; j < x.size(); ++j){
            temp.push_back(pow(x[i], j));
        }
        X.push_back(temp);
        Y.push_back({y[i]});
    }

    std::vector<std::vector<double>> coef = MMULT(INVERSE(X), Y);

    std::vector<double> Fx, Fy;
    int steps = 50;
    double x0 = 1;
    double x1 = 10;
    double dx = (x1 - x0)/((double) steps - 1);

    for(int i = 0; i < steps; ++i){
        double summation = 0;
        double uxp = x0 + i*dx;
        for(int j = 0; j < coef.size(); ++j){
            summation += coef[j][0]*pow(uxp, j);
        }
        Fx.push_back(uxp);
        Fy.push_back(summation);
    }

    plt::scatter2D(ax, x, y, "red");
    plt::plot2D(ax, Fx, Fy, "blue");

    plt::show();

    return 0;
}