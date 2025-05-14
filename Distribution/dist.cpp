#include <iostream>
#include <string>
#include <thread>
#include <chrono>
#include <time.h>
#include <vector>
#include <map>
#include <curl/curl.h>
#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Python.h>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace boost::property_tree;

// Importing stock data from Polygon.io
std::string address(std::string ticker){
    std::string key = "";
    std::string url = "https://api.polygon.io/v2/aggs/ticker/" + ticker + "/range/1/day/2024-01-09/2024-09-29?adjusted=true&sort=asc&limit=300&apiKey=" + key;
    return url;
}

size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::string* response) {
    size_t totalSize = size * nmemb;
    response->append((char*)contents, totalSize);
    return totalSize;
}

std::string RequestData(const std::string& url) {
    CURL* curl;
    CURLcode res;
    std::string response;

    curl = curl_easy_init(); // Initialize CURL
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str()); // Set the URL
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback); // Set the callback function
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response); // Pass the response string to the callback
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L); // Follow redirects if necessary

        // Perform the request
        res = curl_easy_perform(curl);
        
        // Check for errors
        if(res != CURLE_OK) {
            std::cerr << "CURL request failed: " << curl_easy_strerror(res) << std::endl;
        }

        // Cleanup
        curl_easy_cleanup(curl);
    }

    return response; // Return the response
}

void Sleep(int wait_time){
    std::this_thread::sleep_for(std::chrono::seconds(wait_time));
}

std::map<std::string, std::vector<double>> ImportHistoricalData(std::vector<std::string> tickers, int wait_time){
    std::map<std::string, std::vector<double>> result;
    std::cout << "Stock data is loading" << std::endl;
    for(auto & stock : tickers){
        Sleep(wait_time);
        std::cout << stock << " is imported" << std::endl;
        std::string response = RequestData(address(stock));
        ptree data;
        std::stringstream ss(response);
        read_json(ss, data);
        for(ptree::const_iterator it = data.begin(); it != data.end(); ++it){
            if(it->first == "results"){
                for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                    for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                        if(kt->first == "l"){
                            result[stock].push_back(atof(kt->second.get_value<std::string>().c_str()));
                        }
                    }
                }
            }
        }
    }

    return result;
}

void ComputeData(std::vector<double> close, std::string ticker, std::map<std::string, std::map<std::string, std::vector<double>>> & modify){
    auto average = [](std::vector<double> x){
        double total = 0;
        for(auto & i : x){
            total += i;
        }
        total /= x.size();
        return total;
    };

    auto volatility = [&](std::vector<double> x){
        double mu = average(x);
        double total = 0;
        for(auto & i : x){
            total += pow(i - mu, 2);
        }
        total /= ((double) x.size() - 1);
        return pow(total, 0.5);
    };

    auto dWT = [](){
        int num = 10;
        double dw = (rand() % (2*num + 1)) - num;
        return dw/100.0;
    };

    auto histogram = [](std::vector<double> returns, int bins){
        std::map<std::string, std::vector<double>> res;
        std::sort(returns.begin(), returns.end());
        double m0 = returns[0];
        double m1 = returns[returns.size() - 1];
        double dm = (m1 - m0)/((double) bins);
        for(int i = 0; i < bins; ++i){
            double a = m0 + i*dm;
            double b = m0 + (i+1)*dm;
            int count = 0;
            for(auto & r : returns){
                if(i == bins - 1){
                    if(r >= a && r <= b){
                        count += 1;
                    }
                } else {
                    if(r >= a && r < b){
                        count += 1;
                    }
                }
            }
            double mid = 0.5*(a + b);
            res["x"].push_back(mid);
            res["y"].push_back(count);
        }
        return res;
    };

    std::vector<double> ror, ror_predict, stock_paths;
    for(int i = 1; i < close.size(); ++i){
        ror.push_back(close[i]/close[i-1] - 1.0);
    }
    double S = close[close.size() - 1];
    double mu = average(ror);
    double t = 1.0/12.0;
    double v = volatility(ror);
    int N = 1000;
    int P = 100;
    double dt = t / (double) N;

    // mu*S*dt + v*S*dWT()
    srand(time(NULL));

    for(int p = 0; p < P; ++p){
        double S0 = S;
        for(int t = 0; t < N; ++t){
            S0 += mu*S0*dt + v*S0*dWT();
        }
        stock_paths.push_back(S0);
    }

    for(int i = 1; i < stock_paths.size(); ++i){
        ror_predict.push_back(stock_paths[i]/stock_paths[i-1] - 1.0);
    } 

    std::map<std::string, std::vector<double>> RHist = histogram(ror, 30);
    std::map<std::string, std::vector<double>> RPred = histogram(ror_predict, 30);

    modify[ticker]["xhist"] = RHist["x"];
    modify[ticker]["yhist"] = RHist["y"];
    modify[ticker]["xpred"] = RPred["x"];
    modify[ticker]["ypred"] = RPred["y"];

}


int main()
{
    std::vector<std::string> tickers = {"MSFT","AAPL","NVDA","AMZN","IBM","ORCL"};
    std::vector<PyObject*> plots;
    std::vector<int> pnum = {231, 232, 233, 234, 235, 236};
    for(auto & number : pnum){
        plots.push_back(plt::chart2D(number));
    }

    int sleep_for_time = 10;

    std::map<std::string, std::vector<double>> close = ImportHistoricalData(tickers, sleep_for_time);
    std::map<std::string, std::map<std::string, std::vector<double>>> modify;
    
    std::vector<std::thread> items;
    for(auto & ticker : tickers){
        items.emplace_back(ComputeData, close[ticker], ticker, std::ref(modify));
    }

    for(auto & plane : items){
        plane.join();
    }

    for(int i = 0; i < plots.size(); ++i){
        PyObject * ax = plots[i];
        std::string tick = tickers[i];
        plt::PlotTitle(ax, tick);
        plt::plot2D(ax, modify[tick]["xhist"], modify[tick]["yhist"], "red");
        plt::plot2D(ax, modify[tick]["xpred"], modify[tick]["ypred"], "limegreen");
    }

    plt::show();

    return 0;
}