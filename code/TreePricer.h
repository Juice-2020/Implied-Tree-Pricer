#ifndef _TREE_PRICER
#define _TREE_PRICER

#include <vector>
#include "TreeProduct.h"
#include <string>
#include <cmath>
#include "Spline.h"
#include "FunctionNeed.h"
namespace ublas = boost::numeric::ublas;

class BinomialTreePricer
{
 public:
  BinomialTreePricer(int N) {
    nTimeSteps = N;
    states.resize(N+1);
  }
  double Price(const Market& mkt, const TreeProduct& trade)
  {
    // model setup
    double T = trade.GetExpiry() - mkt.asof;
    double dt = T / nTimeSteps;
    ModelSetup(mkt.stockPrice, mkt.volatility, mkt.interestRate, dt);
    // initialize
    for (int i = 0; i <= nTimeSteps; i++) {
      states[i] = trade.Payoff( GetSpot(nTimeSteps, i) );
    }    
    // price by backward induction
    for (int k = nTimeSteps-1; k >= 0; k--)
      for (int i = 0; i <= k; i++) {
	// calculate continuation value
	double df = exp(-mkt.interestRate*dt);	  
	double continuation = df * (states[i]*GetProbUp() + states[i+1]*GetProbDown());
	// calculate the option value at node(k, i)
	states[i] = trade.ValueAtNode( GetSpot(k, i), dt*k, continuation);
      }
    return states[0];
  }
  
 protected:
  virtual void ModelSetup(double S0, double sigma, double rate, double dt) = 0;
  virtual double GetSpot(int ti, int si) const = 0;
  virtual double GetProbUp() const = 0;
  virtual double GetProbDown() const = 0;
  
 private:
  int nTimeSteps;
  std::vector<double> states;
};

class CRRBinomialTreePricer : public BinomialTreePricer
{
 public:
  CRRBinomialTreePricer(int N) : BinomialTreePricer(N){}
 protected:
  virtual void ModelSetup(double S0, double sigma, double rate, double dt)
  {
    double b = std::exp((2*rate+sigma*sigma)*dt)+1;
    u = (b + std::sqrt(b*b - 4*std::exp(2*rate*dt))) / 2 / std::exp(rate*dt);
    p = (std::exp(rate*dt) -  1/u) / (u - 1/u);
    currentSpot = S0;
  }
  virtual double GetSpot(int ti, int si) const
  {
    return currentSpot * std::pow(u, ti-2*si); 
  }
  virtual double GetProbUp() const {return p;}
  virtual double GetProbDown() const {return 1-p;}
  
 private:
  double u; // up multiplicative
  double p; // probability for up state
  double currentSpot; // current market spot price
};

class JRRNBinomialTreePricer : public BinomialTreePricer
{
 public:
  JRRNBinomialTreePricer(int N) : BinomialTreePricer(N){}
 protected:
  virtual void ModelSetup(double S0, double sigma, double rate, double dt)
  {
    u = std::exp( (rate - sigma*sigma/2) * dt + sigma * std::sqrt(dt) );
    d = std::exp( (rate - sigma*sigma/2) * dt - sigma * std::sqrt(dt) );
    p = (std::exp(rate*dt) -  d) / (u - d);
    currentSpot = S0;
  }
  virtual double GetSpot(int ti, int si) const
  {
    return currentSpot * std::pow(u, ti-si) * std::pow(d, si);
  }
  virtual double GetProbUp() const {return p;}
  virtual double GetProbDown() const {return 1-p;}

 private:
  double u; // up multiplicative
  double d; // down multiplicative
  double p; // probability for up state
  double currentSpot; // current market spot price
}; 


// plain binimoal tree pricer
double binomialTree ( const Market& market, const TreeProduct& trade, int N)
{
  // set up the tree model and parameters
  double T = trade.GetExpiry() - market.asof;
  double sigma = market.volatility;
  double rate = market.interestRate;
  std::vector<double> states(N+1);
  double dt = T / N;
  double b = std::exp((2*rate+sigma*sigma)*dt)+1;
  double u = (b + std::sqrt(b*b - 4*std::exp(2*rate*dt))) / 2 / std::exp(rate*dt);
  double p = (std::exp(rate*dt) -  1/u) / (u - 1/u);
  // initialize the final states, apply payoff directly
  for (int i = 0; i <= N; i++) {
    double S = market.stockPrice * std::pow(u, N-2*i);
    states[i] = trade.Payoff(S);
  }
  for (int k = N-1; k >= 0; k--)
    for (int i = 0; i <= k; i++) {
      // calculate continuation value
      double df = exp(-rate*dt);
      double continuation = df * (states[i]*p + states[i+1]*(1-p));
      // calculate the option value at node(k, i)
      double S = market.stockPrice * std::pow(u, k-2*i);
      states[i] = trade.ValueAtNode(S, dt*k, continuation);
    }
  return states[0];
}



double impliedTree_pricer(const Market& market, const TreeProduct& trade, int N, const ImpVol& impvol,double q = 0)
{

    // set up the tree model and parameters
    double T = trade.GetExpiry() - market.asof;
    double rate = market.interestRate;
    double S0 = market.stockPrice;
    std::vector<double> states(N + 1);
    double dt = T / N;
    // initialize the final states, apply payoff directly
    Result set = ImpliedTree(impvol, T, S0, rate, q, N);
    for (int i = 0; i <= N; i++) {
        states[i] = trade.Payoff(set.stock(N, i));
    }
    for (int k = N - 1; k >= 0; k--)
        for (int i = 0; i <= k; i++) {
            // calculate continuation value
            double df = exp(-rate * dt);
            double p = set.p(k+1, i);
            //cout << "p:" << p << "at" << k << "at" << i << endl;
            double continuation = df * (states[i] * p + states[i + 1] * (1 - p));
            //cout << "zhe:" << states[i] * p + states[i + 1] * (1 - p) << endl;
            //cout << "rate" << df;
            //cout << "value:" << continuation << endl;
            // calculate the option value at node(k, i)
            double S = set.stock(k, i);
            states[i] = trade.ValueAtNode(S, dt * k, continuation);
            //cout << "node:" << states[i] << endl;
        }
    return states[0];
}

#endif
