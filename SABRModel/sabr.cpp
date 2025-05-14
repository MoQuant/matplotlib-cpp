#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <algorithm>
#include <curl/curl.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "matplotlibcpp.h"
#include <Python.h>

namespace plt = matplotlibcpp;

using namespace boost::property_tree;

std::string fmp_address(std::string ticker){
    std::string url = "https://financialmodelingprep.com";
    std::string key = "";
    std::string endpoint = "/api/v3/historical-price-full/" + ticker + "?apikey=" + key;
    return url + endpoint;
}

// Callback function to handle the data received from the GET request
size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* s) {
    size_t newLength = size * nmemb;
    try {
        s->append((char*)contents, newLength);
    } catch (std::bad_alloc& e) {
        // Handle memory problem if needed
        return 0;
    }
    return newLength;
}

// Function to perform a GET request
std::string Request(const std::string& url) {
    CURL* curl;
    CURLcode res;
    std::string readBuffer;

    curl = curl_easy_init();  // Initialize cURL
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());  // Set the URL
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);  // Set the callback function
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &readBuffer);  // Set the buffer to store the response
        res = curl_easy_perform(curl);  // Perform the request

        if(res != CURLE_OK) {
            std::cerr << "cURL error: " << curl_easy_strerror(res) << std::endl;
        }

        curl_easy_cleanup(curl);  // Clean up cURL
    }

    return readBuffer;
}

std::vector<double> PullStockData(std::string ticker){
    std::vector<double> close;
    std::string response = Request(fmp_address(ticker));
    std::stringstream ss(response);
    ptree dataset;
    read_json(ss, dataset);
    for(ptree::const_iterator it = dataset.begin(); it != dataset.end(); ++it){
        if(it->first == "historical"){
            for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                    if(kt->first == "adjClose"){
                        close.push_back(kt->second.get_value<double>());
                    }
                }
            }
        }
    }
    std::reverse(close.begin(), close.end());
    return close;
}

double stockPrice(std::vector<double> close){
    return close[close.size() - 1];
}

std::map<std::string, double> Parameters(std::vector<double> close){
    auto mean = [](std::vector<double> x){
        double average = 0;
        for(auto & i : x){
            average += i;
        }
        return average / (double) x.size();
    };
    
    auto stdev = [&](std::vector<double> x){
        double mu = mean(x);
        double volatility = 0;
        for(int i = 0; i < x.size(); ++i){
            volatility += pow(x[i] - mu, 2);
        }
        return pow(volatility / ((double) x.size() - 1), 0.5);
    };
    std::map<std::string, double> result;
    std::vector<double> ror;
    for(int i = 1; i < close.size(); ++i){
        ror.push_back(close[i]/close[i-1] - 1.0);
    }
    double ivol = stdev(ror);
    double window = 100;
    std::vector<double> store_vol;
    for(int i = window; i < ror.size(); ++i){
        std::vector<double> hold = {ror.begin()+i-window, ror.begin()+i};
        store_vol.push_back(stdev(hold));
    }
    double vvol = stdev(store_vol);
    result["iv"] = ivol;
    result["sv"] = vvol;
    return result;
}

std::map<std::string, std::vector<std::vector<double>>> GRID(std::vector<double> x, std::vector<double> y, std::vector<double> T){
    std::map<std::string, std::vector<std::vector<double>>> result;
    for(int i = 0; i < x.size(); ++i){
        std::vector<double> tempx;
        std::vector<double> tempy;
        std::vector<double> tempz;
        for(int j = 0; j < y.size(); ++j){
            tempx.push_back(x[i]);
            tempy.push_back(y[j]);
            tempz.push_back(T[j]);
        }
        result["Strikes"].push_back(tempx);
        result["Forward"].push_back(tempy);
        result["Expiry"].push_back(tempz);
    }
    return result;
}

double ImpliedVol(double iv, double vvol, double rho, double Ft, double K, double beta){
    double z = (vvol/iv)*pow(Ft*K, (1 - beta)/2.0)*log(Ft/K);
    double xz = log((pow(1 - 2.0*rho*z + pow(z, 2), 0.5) + z - rho)/(1 - rho));
    double adj = (pow(1 - beta, 2)/24.0)*(pow(iv, 2)/pow(Ft*K, 1-beta))+(rho*beta*iv*vvol)/(4.0*pow(Ft*K, (1-beta)/2.0))+(pow(vvol, 2)*(2 - 3*pow(rho, 2)))/24.0;
    return (iv/pow(Ft*K, (1 - beta)/2.0))*(z/xz)*(1+adj);
}

std::vector<std::vector<double>> VolSurface(std::map<std::string, std::vector<std::vector<double>>> xy, double iv, double vvol, double rho, double beta){
    std::vector<std::vector<double>> vs;
    std::vector<double> ts;
    for(int i = 0; i < xy["Strikes"].size(); ++i){
        ts.clear();
        for(int j = 0; j < xy["Forward"].size(); ++j){
            ts.push_back(ImpliedVol(iv, vvol, rho, xy["Forward"][i][j], xy["Strikes"][i][j], beta));
        }
        vs.push_back(ts);
    }
    return vs;
}

std::vector<double> linspace(double a, double b, int n){
    std::vector<double> result;
    double dx = (b - a)/(n - 1);
    for(int i = 0; i < n; ++i){
        result.push_back(a + i*dx);
    }
    return result;
}

int main()
{
    std::vector<std::string> tickers = {"MSFT","AAPL","AMZN","IBM","NVDA","GOOGL"};
    std::vector<PyObject*> bx;
    for(auto & num : {231, 232, 233, 234, 235, 236}){
        bx.push_back(plt::chart(num));
    }

    double rho = -0.9;
    double r = 0.044;

    for(int k = 0; k < tickers.size(); ++k){
        std::string ticker = tickers[k];
        PyObject * ax = bx[k];
        std::vector<double> close = PullStockData(ticker);
        double S = stockPrice(close);

        std::map<std::string, double> VOL = Parameters(close);

        std::vector<double> T = linspace(7/365.0, 2.0, 100);
        std::vector<double> Fr;
        for(int t = 0; t < T.size(); ++t){
            Fr.push_back(S*exp(r*T[t]));
        }
        
        std::vector<double> K = linspace(0.5*S, 1.5*S, 100);
        
        double beta = 0.3;
        std::map<std::string, std::vector<std::vector<double>>> grid = GRID(K, Fr, T);
        std::vector<std::vector<double>> ZVol = VolSurface(grid, VOL["iv"], VOL["sv"], rho, beta);

        plt::PlotTitle(ax, ticker);
        plt::surface3DMap(ax, grid["Strikes"], grid["Expiry"], ZVol, "jet", 0.8);
        plt::Chart3DAxesNames(ax, "Strike Price", "Expiry", "Implied Volatility");
    }

    plt::show();

    return 0;
}

/* PART ONE
int main()
{
    std::string ticker = "MSFT";
    double rho = -0.5;
    double r = 0.044;

    std::vector<double> close = PullStockData(ticker);
    double S = stockPrice(close);

    std::map<std::string, double> VOL = Parameters(close);

    std::vector<double> T = {7/365.0, 14/365.0, 30/365.0, 60/365.0, 90/365.0, 1.0, 1.5, 2.0};
    std::vector<double> Fr;
    for(int t = 0; t < T.size(); ++t){
        Fr.push_back(S*exp(r*T[t]));
    }
    
    std::vector<double> K = {390, 400, 410, 420, 430, 440, 450, 460};
    PyObject * ax = plt::chart(111);

    for(double beta = 0.0; beta <= 1.0; beta += 0.01){
        std::map<std::string, std::vector<std::vector<double>>> grid = GRID(K, Fr, T);
        std::vector<std::vector<double>> ZVol = VolSurface(grid, VOL["iv"], VOL["sv"], rho, beta);

        plt::Clear3DChart(ax);
        plt::surface3DMap(ax, grid["Strikes"], grid["Expiry"], ZVol, "jet", 0.8);
        plt::Chart3DAxesNames(ax, "Strike Price", "Expiry", "Implied Volatility");
        plt::pause(0.001);
    }
    plt::show();

    return 0;
}
*/