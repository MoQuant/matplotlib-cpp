// premium += notional*cds*dt/pow(1+r,i+1)
// protection += notional*lgd*p[i]/pow(1 + r, i+1)
// x = lgd y = default

#include <iostream>
#include <vector>
#include <Python.h>
#include <math.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

double premium_leg(double notional, double cds, double dt, double r, int n){
    double pv = 0;
    for(int i = 0; i < n; ++i){
        pv += notional*cds*dt/pow(1 + r, i + 1);
    }
    return pv;
}

double dpremium_leg(double notional, double dt, double r, int n){
    double pv = 0;
    for(int i = 0; i < n; ++i){
        pv += notional*dt/pow(1 + r, i + 1);
    }
    return pv;
}

double protection_leg(double notional, double lgd, double r, std::vector<double> pb){
    double pv = 0;
    for(int i = 0; i < pb.size(); ++i){
        pv += notional*lgd*pb[i]/pow(1 + r, i + 1);
    }
    return pv;
}

double calculate_cds(double notional, double lgd, double dt, double r, int n, std::vector<double> pb){
    double cds_0 = 0.01, cds_1 = 0.99;
    while(true){
        cds_1 = cds_0 - (premium_leg(notional, cds_0, dt, r, n) - protection_leg(notional, lgd, r, pb))/dpremium_leg(notional, dt, r, n);
        if(fabs(cds_1 - cds_0) <= 0.00001){
            break;
        }
        cds_0 = cds_1;
    }
    return cds_1;
}

std::vector<std::vector<double>> transpose(std::vector<std::vector<double>> x){
    std::vector<std::vector<double>> z;
    std::vector<double> temp;
    for(int i = 0; i < x[0].size(); ++i){
        temp.clear();
        for(int j = 0; j < x.size(); ++j){
            temp.push_back(x[j][i]);
        }
        z.push_back(temp);
    }
    return z;
}

std::vector<double> Reverse(std::vector<double> x){
    std::vector<double> res;
    double product = 1.0;
    res.push_back(1.0);
    for(int i = 0; i < x.size(); ++i){
        product *= (1 - x[i]);
        res.push_back(product);
    }
    return res;
}

int main()
{
    std::vector<double> the_time = {30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0};
    double notional = 1000000;
    PyObject * ax = plt::chart(111);

    for(auto & ts : the_time){
        double t = ts/365.0;
        int n = 50;
        double dt = t / (double) n;
        double rf = pow(1 + 0.0548, 1/30) - 1;

        std::vector<double> default_rates = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09};
        std::vector<double> recovery_rates = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        std::reverse(recovery_rates.begin(), recovery_rates.end());
        
        default_rates = Reverse(default_rates);

        std::vector<std::vector<double>> x, y, z;
        std::vector<double> temp;
        
        for(int i = 0; i < default_rates.size(); ++i){
            x.push_back(recovery_rates);
            y.push_back(default_rates);
        }
        y = transpose(y);

        for(int i = 0; i < x.size(); ++i){
            temp.clear();
            for(int j = 0; j < x[0].size(); ++j){
                double lgd = x[i][j];
                std::vector<double> pb = y[i];
                double calculation = calculate_cds(notional, lgd, dt, rf, n, pb);
                temp.push_back(calculation);   
            }
            z.push_back(temp);
        }

        
        plt::Chart3DAxesNames(ax, "LGD", "Probability", "CDS Spread");

        plt::surface3DMap(ax, x, y, z, "jet", 1.0);

        plt::pause(1);
    }

    plt::show();


    return 0;
}