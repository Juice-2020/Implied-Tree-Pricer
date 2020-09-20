#ifndef _TREE_PRODUCT_H
#define _TREE_PRODUCT_H
#include "Date.h"

class TreeProduct
{
 public:
  virtual const Date& GetExpiry() const = 0;
  virtual double Payoff(double stockPrice) const = 0;
  virtual double ValueAtNode(double stockPrice,
			     double t,
			     double continuationValue) const = 0;
};
#endif
