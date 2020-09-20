#include "Payoff.h"
#include <iostream>
#include <cmath> // for std::exp()
#include <cassert> // for assertion on inputs
#include "bsformula.h"
double CRR(OptType optType, double K, double T, double S_0, double sigma, double rate, double q, int N)
{
    int nTimeSteps = N;
    double dt = T / N;
    double b = std::exp((2 * rate + sigma * sigma) * dt) + 1;
    double u = (b + std::sqrt(b * b - 4 * std::exp(2 * rate * dt))) / 2 / std::exp(rate * dt);
    double p = (std::exp(rate * dt) - 1 / u) / (u - 1 / u);
    std::vector<double> states(N + 1);
    switch (optType)
    {
    case C:
        for (int i = 0; i <= N; i++) {
            double S = S_0 * std::pow(u, N - 2 * i);
            states[i] = PAYOFF::VanillaOption(Call, K, S);
        }
        for (int k = N - 1; k >= 0; k--)
            for (int i = 0; i <= k; i++) {
                // calculate continuation value
                double df = exp(-rate * dt);
                double continuation = df * (states[i] * p + states[i + 1] * (1 - p));
                // calculate the option value at node(k, i)
                double S = S_0 * std::pow(u, k - 2 * i);
                states[i] = continuation;
            }
        return states[0];
        break;
    case P:
        for (int i = 0; i <= N; i++) {
            double S = S_0 * std::pow(u, N - 2 * i);
            states[i] = PAYOFF::VanillaOption(Put, K, S);
        }
        for (int k = N - 1; k >= 0; k--)
            for (int i = 0; i <= k; i++) {
                // calculate continuation value
                double df = exp(-rate * dt);
                double continuation = df * (states[i] * p + states[i + 1] * (1 - p));
                // calculate the option value at node(k, i)
                double S = S_0 * std::pow(u, k - 2 * i);
                states[i] = continuation;
            }
        return states[0];
        break;
    default:
        throw "unsupported optionType";
    }

    // std::cout << "The--"<< optType << " --option price is " << V_0 << std::endl;
}
