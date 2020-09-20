#ifndef _EUROPEAN_TRADE
#define _EUROPEAN_TRADE

#include <cassert> 

#include "TreeProduct.h"
#include "Payoff.h"
#include "Types.h"

// still a abstract class
class EuropeanTrade : public TreeProduct
{
 public:
  virtual double ValueAtNode(double S, double t, double continuation) const
  { return continuation; }
};

class EuropeanOption : public EuropeanTrade
{
 public:
  virtual double Payoff(double S) const
  {
    return PAYOFF::VanillaOption(optType, strike, S);
  }
  virtual const Date& GetExpiry() const
  { return expiryDate; }
  EuropeanOption(OptionType _optType, double _strike, const Date& _expiry)
    : optType(_optType), strike(_strike), expiryDate(_expiry) {}
 private:
  OptionType optType;
  double strike;
  Date expiryDate;  
};

class EuroCallSpread : public EuropeanTrade
{
 public:
  virtual double Payoff(double S) const
  {
    return PAYOFF::CallSpread(strike1, strike2, S);
  }
  virtual const Date& GetExpiry() const
  { return expiryDate; }
  EuroCallSpread(double _k1, double _k2, const Date& _expiry)
   : strike1(_k1), strike2(_k2), expiryDate(_expiry)
  {
    assert(_k1 < _k2);
  };
 private:
  double strike1;
  double strike2;
  Date expiryDate;
};


#endif
