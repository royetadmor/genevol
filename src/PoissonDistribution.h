#ifndef POISSON_DISTRIBUTION_H
#define POISSON_DISTRIBUTION_H

using namespace std;

class PoissonDistribution
{
public:
    /**
     * Construct a truncated Poisson distribution.
     */
    PoissonDistribution(double lambda, int maxK);

    /**
     * Return the normalized probability P(X = k).
     * Returns 0 for invalid k.
     */
    double probability(int k) const;

    /**
     * Return log P(X = k).
     */
    double logProbability(int k) const;
    double lambda() const noexcept { return lambda_; };
    int maxK() const noexcept { return maxK_; };

private:
    /// Log of the unnormalized Poisson PMF
    double logPoissonPMF(int k) const;

    /// Compute log normalization constant for truncated case
    void computeNormalization();

private:
    double lambda_;
    int maxK_;  
    double logZ_;
};

#endif // POISSON_DISTRIBUTION_H
