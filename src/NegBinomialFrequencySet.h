// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef NEGBINOMIAL_FREQUENCY_SET_H
#define NEGBINOMIAL_FREQUENCY_SET_H

#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Numeric/Parameter.h>

namespace bpp
{

/**
 * @brief Root frequency set driven by a truncated negative binomial distribution.
 *
 * The frequency of state k is:
 *   P(X=k) = Γ(k+r) / (k! · Γ(r)) · p^r · (1-p)^k  (normalized over [0, maxK))
 *   with p = r / (r + μ)
 *
 * The dispersion parameter r is registered as a bpp parameter and is intended
 * to be numerically optimized. The empirical mean μ is fixed from the data.
 */
class NegBinomialFrequencySet :
  public AbstractFrequencySet
{
private:
  double mu_;  // fixed empirical mean; p = r / (r + mu)

public:
  /**
   * @param stateMap  State map from the substitution model.
   * @param r         Initial dispersion parameter (> 0).
   * @param mu        Empirical mean of gene family sizes (> 0); fixed.
   */
  NegBinomialFrequencySet(
    std::shared_ptr<const StateMapInterface> stateMap,
    double r,
    double mu);

  NegBinomialFrequencySet(const NegBinomialFrequencySet&) = default;

  NegBinomialFrequencySet* clone() const override
  {
    return new NegBinomialFrequencySet(*this);
  }

public:
  std::string getName() const override { return "NegBinomial"; }

  void setFrequencies(const std::vector<double>& frequencies) override;

  double getR()  const { return getParameterValue("r"); }
  double getMu() const { return mu_; }

  void printFrequencies() const;

protected:
  void fireParameterChanged(const ParameterList& parameters) override;

private:
  void updateFrequencies_();
};

} // namespace bpp

#endif // NEGBINOMIAL_FREQUENCY_SET_H
