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

// Fetches the stock price data from Financial Modeling Prep
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

// Matrix Multiplication Function
std::vector<std::vector<double>> MMULT(std::vector<std::vector<double>> x,
                                       std::vector<std::vector<double>> y)
{
    std::vector<std::vector<double>> result;
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
        result.push_back(temp);
    }
    return result;
}

// Matrix Transpose Function
std::vector<std::vector<double>> TRANSPOSE(std::vector<std::vector<double>> z)
{
    std::vector<std::vector<double>> X;
    std::vector<double> temp;
    for(int i = 0; i < z[0].size(); ++i){
        temp.clear();
        for(int j = 0; j < z.size(); ++j){
            temp.push_back(z[j][i]);
        }
        X.push_back(temp);
    }
    return X;
}

// Inverse Matrix Function using Gaussian Elimination
std::vector<std::vector<double>> INVERSE(std::vector<std::vector<double>> x)
{
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

    double A, B;

    for(int i = 1; i < n; ++i){
        for(int j = 0; j < i; ++j){
            A = x[i][j];
            B = x[j][j];
            for(int k = 0; k < n; ++k){
                x[i][k] = x[i][k] - (A/B)*x[j][k];
                I[i][k] = I[i][k] - (A/B)*I[j][k];
            }
        }
    }

    for(int i = 1; i < n; ++i){
        for(int j = 0; j < i; ++j){
            A = x[j][i];
            B = x[i][i];
            for(int k = 0; k < n; ++k){
                x[j][k] = x[j][k] - (A/B)*x[i][k];
                I[j][k] = I[j][k] - (A/B)*I[i][k];
            }
        }
    }
    
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            I[i][j] = I[i][j] / x[i][i];
        }
    }

    return I;
}

// Multiplies a matrix by a coeffecient
std::vector<std::vector<double>> FACTOR(double a, std::vector<std::vector<double>> x)
{
    for(int i = 0; i < x.size(); ++i){
        for(int j = 0; j < x[0].size(); ++j){
            x[i][j] *= a;
        }
    }
    return x;
}

// Adds or Subtracts a matrix from a matrix with sign = -1 or 1
std::vector<std::vector<double>> ADDSUB(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b, double sign)
{
    for(int i = 0; i < a.size(); ++i){
        for(int j = 0; j < a[0].size(); ++j){
            a[i][j] += sign*b[i][j];
        }
    }
    return a;
}

// Calculates the rate of return matrix of a given matrix of stock prices
std::vector<std::vector<double>> RateOfReturn(std::vector<std::vector<double>> x)
{
    std::vector<std::vector<double>> y;
    std::vector<double> temp;
    for(int i = 1; i < x.size(); ++i){
        temp.clear();
        for(int j = 0; j < x[0].size(); ++j){
            temp.push_back(x[i][j]/x[i-1][j] - 1.0);
        }
        y.push_back(temp);
    }
    return y;
}

// Fetches the stock price data and stores them into a map with the key being the stock ticker and the value being the vector of close prices
std::map<std::string, std::vector<double>> Cyclone(std::vector<std::string> tickA, std::vector<std::string> tickB){
    std::map<std::string, std::vector<double>> prices;

    // Fetches stocks
    for(auto & ticker : tickA){
        // Fetch the stock data
        std::string resp = Request(fmp_address(ticker));
        std::stringstream ss(resp);
        ptree df;
        // Parse string to JSON using Boost
        read_json(ss, df);
        for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
            if(it->first == "historical"){
                for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                    for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                        if(kt->first == "adjClose"){
                            // Store adjusted close prices in map vector based on stock ticker
                            prices[ticker].push_back(atof(kt->second.get_value<std::string>().c_str()));
                        }
                    }
                }
            }
        }
        // Reverse prices so oldest price is first and newest price is last
        std::reverse(prices[ticker].begin(), prices[ticker].end());
    }

    // Fetches hedging instruments
    for(auto & ticker : tickB){
        // Fetches data from Financial Modeling Prep
        std::string resp = Request(fmp_address(ticker));
        std::stringstream ss(resp);
        ptree df;
        // Parses response as JSON with Boost
        read_json(ss, df);
        for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
            if(it->first == "historical"){
                for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                    for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                        if(kt->first == "adjClose"){
                            // Stores adjusted close prices in map  
                            prices[ticker].push_back(atof(kt->second.get_value<std::string>().c_str()));
                        }
                    }
                }
            }
        }
        // Oldest prices first, newest prices last
        std::reverse(prices[ticker].begin(), prices[ticker].end());
    }
    return prices;
}

// Calculates the Minimum Variance Portfolio
std::vector<double> MinVariancePortfolio(std::vector<std::vector<double>> ror, int lookback)
{
    // Generates a set of weights for each inputted portfolio returns snapshot
    auto solver = [](std::vector<std::vector<double>> x){

        // Sets bounds
        int m = x.size();
        int n = x[0].size();
        std::vector<std::vector<double>> mu, cov, temporary;

        // Calculates the mean vector
        for(int i = 0; i < m; ++i){
            mu.push_back({1.0});
        }
        mu = FACTOR(1.0/((double) m), MMULT(TRANSPOSE(mu), x));

        // Subtracts mean from returns matrix
        for(int i = 0; i < m; ++i){
            for(int j = 0; j < n; ++j){
                x[i][j] -= mu[0][j];
            }
        }

        // Calculates the covariance matrix
        cov = FACTOR(1.0/((double) m - 1), MMULT(TRANSPOSE(x), x));

        // Multiplies covariance matrix by 2.0 for the MinVariance matrix
        cov = FACTOR(2.0, cov);
        std::vector<double> ones, weights;
        mu.clear();

        // Adds 1.0 to each column in the matrix
        for(int i = 0; i < n; ++i){
            cov[i].push_back(1.0);
            ones.push_back(1.0);
            mu.push_back({0.0}); // Necessary for matrix multiplication
        }
        // Adds 1.0 to each row in the matrix
        ones.push_back(0.0);
        cov.push_back(ones);
        mu.push_back({1.0});

        // Computes the weights of the portfolio by taking the inverse matrix and multiplying it by zeros
        temporary = MMULT(INVERSE(cov), mu);
        for(int i = 0; i < n; ++i){
            weights.push_back(temporary[i][0]);
        }
        return weights;
    };

    std::vector<double> result;

    // Takes a rolling snapshot of the rate of returns matrix and inputs it into the solver to generate min-variance weights
    for(int i = lookback; i < ror.size(); ++i){
        // Returns snapshot
        std::vector<std::vector<double>> hold_items = {ror.begin() + (i - lookback), ror.begin() + i};

        // Allocated weights
        std::vector<double> weights = solver(hold_items);

        // Summation to generate weighted portfolio from inputted stocks
        double total = 0;
        for(int j = 0; j < ror[0].size(); ++j){
            total += weights[j]*ror[i][j];
        }
        if(std::isinf(total) || std::isnan(total)){
            total = 0;
        }
        result.push_back(total);
    }

    return result;
}

// Uses the Kalman Filter to calculate a hedging ratio along with its test statistic to test significance
std::vector<double> HedgingRatio(std::vector<double> x, std::vector<double> y){
    std::vector<double> result;

    // Used for building a quick vector off just x[i]
    auto bx = [](double ux){
        std::vector<std::vector<double>> result = {{1.0}, {ux}};
        return result;
    };

    // Calculates the mean of a given vector
    auto mean = [](std::vector<double> q){
        double total = 0;
        for(auto & i : q){
            total += i;
        }
        total /= ((double) q.size());
        return total;
    };
    
    std::vector<std::vector<double>> B, Bp, Pp, Q, P, Yp, K, DK, deltaB;
    double R = 0;

    // Set initial parameters to Kalman Filter inputs
    B = {{0.1}, {0.1}};
    Bp = {{0.1}, {0.1}};
    Pp = {{1.0, 0.0}, {0.0, 1.0}};
    Q = {{1.0, 0.0}, {0.0, 1.0}};
    P = {{1.0, 0.0}, {0.0, 1.0}};

    // Kalman Filter Processing
    for(int i = 0; i < x.size(); ++i){
        // Iterate newest beta
        Bp = B;

        // Add Q and P
        Pp = ADDSUB(Q, P, 1.0);

        // Generate Yhat for current point
        Yp = MMULT(TRANSPOSE(Bp), bx(x[i]));

        // Compute the R parameter which is the average residual
        if(i > 2){
            R = 0;
            for(int t = 0; t < i; ++t){
                R += pow(y[t] - (Bp[0][0] + Bp[1][0]*x[t]), 2);
            }
            R = R / ((double) i - 1);
        }

        // Compute Kalman Gain
        K = MMULT(Pp, bx(x[i]));
        DK = MMULT(TRANSPOSE(K), bx(x[i]));

        // Add R to and divide Kalman Gain
        DK[0][0] += R;
        for(int t = 0; t < K.size(); ++t){
            K[t][0] /= DK[0][0];
        }
        // Add beta to error multiplied by Kalman gain
        B = ADDSUB(Bp, FACTOR(y[i] - Yp[0][0], K), 1.0);

        // Update P by subtracting Pp by the multiplication of Kalman gain by X and Pp
        P = ADDSUB(Pp, MMULT(K, MMULT(TRANSPOSE(bx(x[i])), Pp)), -1.0);

        // Update Q by multiplying the error between beta and predicted beta
        deltaB = ADDSUB(B, Bp, -1.0);
        Q = MMULT(deltaB, TRANSPOSE(deltaB));
    }

    // Compute test statistic off the residual sum of squares and the sum squared of x minus its mean
    double rss = 0;
    double bottom = 0;
    double mux = mean(x);
    for(int i = 0; i < x.size(); ++i){
        rss += pow(y[i] - (B[0][0] + B[1][0]*x[i]), 2);
        bottom += pow(x[i] - mux, 2);
    }

    rss /= ((double) x.size() - 2);

    // Generate test statistic and return hedging ratio with it
    double t_stat = B[1][0] / sqrt(rss/bottom);

    result.push_back(B[0][0]);
    result.push_back(B[1][0]);
    result.push_back(t_stat);

    return result;
}


int main()
{
    // Initialize stock and hedging instrument tickers
    std::vector<std::string> port_ticks = {"AAPL","MSFT","NVDA","GOOGL"};
    std::map<std::string, std::string> hedge_items = {
        {"VDC","Vangaurd Consumer Staples ETF"},
        {"IXJ","iShares Global Healthcare ETF"},
        {"SHY","iShares 1-3 Year Treasury Bond ETF"},
        {"GLDM","SPDR Gold MiniShares"}
    };
    std::vector<std::string> hedge_ticks = {"VDC","IXJ","SHY","GLDM"};

    // Fetch all price data
    std::map<std::string, std::vector<double>> prices = Cyclone(port_ticks, hedge_ticks);

    std::vector<std::vector<double>> closePort, closeHedge, rorPort, rorHedge;

    // Split up stocks and hedging instrument prices
    for(auto & tick : port_ticks){
        closePort.push_back(prices[tick]);
    }

    for(auto & tick : hedge_ticks){
        closeHedge.push_back(prices[tick]);
    }

    // Transpose matrices in order to easily compute the rate of returns
    closePort = TRANSPOSE(closePort);
    closeHedge = TRANSPOSE(closeHedge);

    rorPort = RateOfReturn(closePort);
    rorHedge = RateOfReturn(closeHedge);

    // Transpose hedging instruments back so they are easier to loop through for the kalman filter
    rorHedge = TRANSPOSE(rorHedge);

    int lookback = 100;

    // Generate the Minimum Variance Portfolio from the stock data returns
    std::vector<double> PortfolioReturns = MinVariancePortfolio(rorPort, lookback);

    // Calculate regression line x-axis off the portfolio returns min and max
    auto min_p = std::min_element(PortfolioReturns.begin(), PortfolioReturns.end());
    auto max_p = std::max_element(PortfolioReturns.begin(), PortfolioReturns.end());

    double x0 = (double) *min_p;
    double x1 = (double) *max_p;

    int line_length = 100;
    double dX = (x1 - x0)/(line_length - 1);


    std::vector<double> XP;
    for(int i = 0; i < line_length; ++i){
        XP.push_back(x0 + i*dX);
    }

    // Generate plots
    std::vector<PyObject*> plots;
    for(auto & place : {221, 222, 223, 224}){
        plots.push_back(plt::chart2D(place));
    }

    // Generate hedging parameters and construct regression to be plotted
    for(int i = 0; i < hedge_ticks.size(); ++i){
        
        // Adjust hedging data to match portfolio returns dimensions
        rorHedge[i] = {rorHedge[i].begin() + lookback, rorHedge[i].end()};
        std::vector<double> hp = HedgingRatio(PortfolioReturns, rorHedge[i]);
        
        plt::PlotTitle(plots[i], hedge_items[hedge_ticks[i]]);
        std::vector<double> YP;

        for(int t = 0; t < XP.size(); ++t){
            YP.push_back(hp[0] + hp[1]*XP[t]);
        }

        // Plot the points and regression line
        plt::scatter2D(plots[i], PortfolioReturns, rorHedge[i], "red");
        plt::plot2D(plots[i], XP, YP, "blue");

        // Prints out whether each hedging instrument has a significant ratio or not
        std::cout << hedge_ticks[i] << " has a test-statistic of " << hp[2] << std::endl;
    }

    
    plt::show();


    return 0;
}
