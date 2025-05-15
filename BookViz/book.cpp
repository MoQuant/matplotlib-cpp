#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cpprest/ws_client.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <thread>
#include <algorithm>
#include <chrono>
#include "matplotlibcpp.h"
#include <Python.h>

// Declare global namespaces to simplify calling objects
using namespace web;
using namespace web::websockets::client;
using namespace boost::property_tree;
namespace plt = matplotlibcpp;

// Coinbase datafeed class
class datafeed {
    private:
        // Converts a string response to JSON using Boost
        ptree JSON(std::string message){
            std::stringstream ss(message);
            ptree result;
            read_json(ss, result);
            return result;
        }

        // Parses the level2 orderbook for Bitcoin and updates each change in order which the book receives
        void CYCLONE(ptree df, std::map<double, double> & bids, std::map<double, double> & asks){
            bool snapshot = false;
            bool l2update = false;
            for(ptree::const_iterator it = df.begin(); it != df.end(); ++it){
                if(l2update == true && it->first == "changes"){
                    for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                        std::vector<std::string> hold;
                        for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                            hold.push_back(kt->second.get_value<std::string>().c_str());
                        }

                        // Extract the latest price and volume
                        double price = atof(hold[1].c_str());
                        double volume = atof(hold[2].c_str());

                        // Erases part of book if the volume equals zero meaning that the order has been filled or cancelled fully
                        if(hold[0] == "buy"){
                            if(volume == 0){
                                bids.erase(price);
                            } else {
                                bids[price] = volume;
                            }
                        } else {
                            if(volume == 0){
                                asks.erase(price);
                            } else {
                                asks[price] = volume;
                            }
                        }
                    }
                }
                // Parses the bid book in the initial snapshot
                if(snapshot == true && it->first == "bids"){
                    for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                        std::vector<double> hold;
                        for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                            hold.push_back(atof(kt->second.get_value<std::string>().c_str()));
                        }
                        bids[hold[0]] = hold[1];
                    }
                }
                // Parses the ask book in the initial snapshot
                if(snapshot == true && it->first == "asks"){
                    for(ptree::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt){
                        std::vector<double> hold;
                        for(ptree::const_iterator kt = jt->second.begin(); kt != jt->second.end(); ++kt){
                            hold.push_back(atof(kt->second.get_value<std::string>().c_str()));
                        }
                        asks[hold[0]] = hold[1];
                    }
                }

                // Activates whether the message is a snapshot or an update
                if(it->first == "type"){
                    if(it->second.get_value<std::string>() == "l2update"){
                        l2update = true;
                    }
                    if(it->second.get_value<std::string>() == "snapshot"){
                        snapshot = true;
                    }
                }
            }
        }

        

    public:

        // Level2 limit orderbook websocket feed using cpprest
        static void Socket(datafeed dx, std::map<double, double> & bids, std::map<double, double> & asks){
            std::string url = "wss://ws-feed.exchange.coinbase.com";
            std::string msg = "{\"type\":\"subscribe\",\"product_ids\":[\"BTC-USD\"],\"channels\":[\"level2_batch\"]}";

            // Connect to client and send the subscription message
            websocket_client client;
            client.connect(url).wait();
            websocket_outgoing_message outmsg;
            outmsg.set_utf8_message(msg);
            client.send(outmsg);

            while(true){
                client.receive().then([](websocket_incoming_message inmsg){
                    return inmsg.extract_string();
                }).then([&](std::string message){
                    // Retreives response from websocket and parses the data as a snapshot of the book, or an update of the book
                    dx.CYCLONE(dx.JSON(message), std::ref(bids), std::ref(asks));
                }).wait();
            }

            client.close().wait();
        
        }
};

// Sets the examination depth to 80 and pulls the orderbooks price and volume for both bids and asks
std::map<std::string, std::vector<double>> Extract(std::map<double, double> bids, std::map<double, double> asks)
{
    // Set depth of book
    int depth = 80;
    std::map<std::string, std::vector<double>> result;

    // Pulls bid orders and takes the cumulative volume summation
    int count = 0;
    double bidvol = 0;
    for(auto it = bids.rbegin(); it != bids.rend(); ++it){
        bidvol += it->second;
        result["bidPrice"].push_back(it->first);
        result["bidSize"].push_back(bidvol);
        count += 1;
        // Breaks loop once depth limit has been reached
        if(count >= depth){
            break;
        }
    }

    // Pulls ask orders and takes the cumulative volume summation
    count = 0;
    double askvol = 0;
    for(auto it = asks.begin(); it != asks.end(); ++it){
        askvol += it->second;
        result["askPrice"].push_back(it->first);
        result["askSize"].push_back(askvol);
        count += 1;
        // Breaks loop once depth limit has been reached
        if(count >= depth){
            break;
        }
    }

    // Bids must be reveresed in order to be plotted as the lowest bid must be at the beginning of the vector
    std::reverse(result["bidPrice"].begin(), result["bidPrice"].end());
    std::reverse(result["bidSize"].begin(), result["bidSize"].end());
    

    return result;
}

// Builds a vector with a single value for n elements
std::vector<double> push_into(double ii, int n){
    std::vector<double> result;
    for(int i = 0; i < n; ++i){
        result.push_back(ii);
    }
    return result;
}

int main()
{
    // Declare two plots, one for bid and one for ask
    PyObject * ax = plt::chart(121);
    PyObject * ay = plt::chart(122);

    std::map<double, double> bids, asks;
    std::map<std::string, std::vector<double>> preprice;

    // Open datafeed thread
    datafeed wsfeed;
    std::thread feed(wsfeed.Socket, wsfeed, std::ref(bids), std::ref(asks));

    // Sleep for 10 seconds to allow the orderbook thread to build a dataset
    std::this_thread::sleep_for(std::chrono::seconds(10));

    std::vector<std::vector<double>> bX, bY, bZ, aX, aY, aZ;

    int ii = 0;
    int jj = 0;
    int limit = 80;

    while(true){
        // Extract the depth orderbook data
        preprice = Extract(bids, asks);

        // Store bid snapshots for 3D plotting
        bX.clear();
        bY.clear();
        bZ.push_back(preprice["bidSize"]);
        for(int u = 0; u < bZ.size(); ++u){
            bY.push_back(preprice["bidPrice"]);
            bX.push_back(push_into(ii, preprice["bidPrice"].size()));
            ii += 1;
        }

        // Store ask snapshots for 3D plotting
        aX.clear();
        aY.clear();
        aZ.push_back(preprice["askSize"]);
        for(int u = 0; u < aZ.size(); ++u){
            aY.push_back(preprice["askPrice"]);
            aX.push_back(push_into(jj, preprice["askPrice"].size()));
            jj += 1;
        }

        // Erase data from the beginning of the vectors if limit has been reached
        if(bZ.size() >= limit){
            bX.erase(bX.begin());
            bY.erase(bY.begin());
            bZ.erase(bZ.begin());
        }

        if(aZ.size() >= limit){
            aX.erase(aX.begin());
            aY.erase(aY.begin());
            aZ.erase(aZ.begin());
        }

        // Clear and plot the two orderbooks in their own plots
        plt::Clear3DChart(ax);
        plt::Clear3DChart(ay);

        plt::surface3D(ax, bX, bY, bZ, "red", 1.0);
        plt::surface3D(ay, aX, aY, aZ, "limegreen", 1.0);

        plt::pause(0.1);

        //std::cout << bX.size() << " " << bX[0].size() << "\t" << bY.size() << " " << bY[0].size() << "\t" << bZ.size() << " " << bZ[0].size() << std::endl;
    }

    plt::show();
    feed.join();
    return 0;
}
