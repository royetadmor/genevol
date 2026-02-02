#include "PoissonDistribution.h"

#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>


using namespace std;

PoissonDistribution::PoissonDistribution(double lambda, int maxK)
    : lambda_(lambda), maxK_(maxK)
{
    if (lambda_ <= 0.0)
        throw std::invalid_argument("Poisson lambda must be > 0");
    if (maxK_ < 0)
        throw std::invalid_argument("maxK must be >= 0");

    computeNormalization();
}

double logSumExp(const std::vector<double>& logVals)
{
    double maxLog = -std::numeric_limits<double>::infinity();
    for (double v : logVals) 
        if (v > maxLog)
            maxLog = v;

    double sum = 0.0;
    for (double v : logVals)
        sum += std::exp(v - maxLog);

    return maxLog + std::log(sum);
}

double PoissonDistribution::logPoissonPMF(int k) const
{
    // log P(X = k) = k * log(lambda) - lambda - log(k!)
    return k * std::log(lambda_) - lambda_ - std::lgamma(static_cast<double>(k) + 1.0);
}

void PoissonDistribution::computeNormalization()
{
    std::vector<double> logVals;
    logVals.reserve(maxK_);

    for (int k = 0; k < maxK_; ++k)
        logVals.push_back(logPoissonPMF(k));

    logZ_ = logSumExp(logVals);
}

double PoissonDistribution::logProbability(int k) const
{
    if (k < 0)
        return -std::numeric_limits<double>::infinity();

    if (maxK_ >= 0)
    {
        if (k >= maxK_)
            return -std::numeric_limits<double>::infinity();

        return logPoissonPMF(k) - logZ_;
    }

    // Untruncated case
    return logPoissonPMF(k);
}

double PoissonDistribution::probability(int k) const
{
    return std::exp(logProbability(k));
}