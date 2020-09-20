//
//  BSformula.h
//  Implied Tree
//
//  Created by 陈家豪 on 2020/4/13.
//  Copyright © 2020 陈家豪. All rights reserved.
//

#ifndef BSformula_h
#define BSformula_h
#include <iostream>
#include <cmath> // for std::exp()
#include <cassert> // for assertion on inputs
enum OptType{C, P};

double cnorm(double x);

double bsformula(OptType optType,double K,double T, double S_0,double sigma,double rate,double q)
{
  double sigmaSqrtT = sigma * std::sqrt(T);
  double d1 = (std::log(S_0 / K) + (rate-q) * T)/sigmaSqrtT + 0.5 * sigmaSqrtT;
  double d2 = d1 - sigmaSqrtT;
    
  double V_0;
  switch (optType)
    {
    case C:
      V_0 = S_0 *exp(-q*T)* cnorm(d1) - K * exp(-rate*T) * cnorm(d2);
      break;
    case P:
      V_0 = K * exp(-rate*T) * cnorm(-d2) - S_0 *exp(-q*T)* cnorm(-d1);
      break;
    default:
      throw "unsupported optionType";
    }

 // std::cout << "The--"<< optType << " --option price is " << V_0 << std::endl;
  return V_0;
}

double cnorm(double x)
{
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;
  int sign = 1;
  if (x < 0)
    sign = -1;
  x = fabs(x)/sqrt(2.0);
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
  return 0.5*(1.0 + sign*y);
}


#endif /* BSformula_h */
