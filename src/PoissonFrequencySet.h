// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef POISSON_FREQUENCY_SET_H
#define POISSON_FREQUENCY_SET_H

#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>

namespace bpp
{

/**
 * @brief Root frequency set driven by a truncated Poisson distribution.
 *
 * The frequency of state k is:
 *   P(X=k) = e^(-lambda) * lambda^k / k!  (normalized over [0, maxK))
 *
 */
class PoissonFrequencySet :
  public AbstractFrequencySet
{
private:
  double lambda_;
  bool optimize_;

public:
  /**
   * @param stateMap  State map from the substitution model.
   * @param lambda    Mean of the Poisson distribution (starting value when optimize=true).
   * @param optimize  If true, lambda is registered as a bpp Parameter and optimized.
   */
  PoissonFrequencySet(
    std::shared_ptr<const StateMapInterface> stateMap,
    double lambda,
    bool optimize = false);

  PoissonFrequencySet(const PoissonFrequencySet&) = default;

  PoissonFrequencySet* clone() const override
  {
    return new PoissonFrequencySet(*this);
  }

public:
  std::string getName() const override { return "Poisson"; }

  void setFrequencies(const std::vector<double>& frequencies) override;

  double getLambda() const { return lambda_; }

  void fireParameterChanged(const ParameterList& parameters) override;

private:
  void updateFrequencies_();
};

}

#endif // POISSON_FREQUENCY_SET_H