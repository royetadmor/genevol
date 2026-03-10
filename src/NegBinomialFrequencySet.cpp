// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "NegBinomialFrequencySet.h"

#include <Bpp/Numeric/Parameter.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace bpp;
using namespace std;

NegBinomialFrequencySet::NegBinomialFrequencySet(
  shared_ptr<const StateMapInterface> stateMap,
  double r,
  double mu)
  :
  AbstractFrequencySet(stateMap, "NegBinomial.", "NegBinomial"),
  mu_(mu)
{
  if (r <= 0.0)
    throw Exception("NegBinomialFrequencySet: r must be > 0.");
  if (mu_ <= 0.0)
    throw Exception("NegBinomialFrequencySet: mu must be > 0.");

  addParameter_(new Parameter("NegBinomial.r", r, Parameter::R_PLUS_STAR));

  updateFrequencies_();
}

void NegBinomialFrequencySet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequencySet::fireParameterChanged(parameters);
  updateFrequencies_();
}

void NegBinomialFrequencySet::printFrequencies() const
{
  const std::vector<double>& freqs = getFrequencies();
  std::cout << "NegBinomial root frequencies (r=" << getR() << ", mu=" << mu_ << ", p=" << getR() / (getR() + mu_) << "):" << std::endl;
  for (size_t k = 0; k < freqs.size(); ++k)
    std::cout << "  P(X=" << k << ") = " << freqs[k] << std::endl;
}

void NegBinomialFrequencySet::setFrequencies(const std::vector<double>& frequencies)
{
  if (frequencies.size() != getNumberOfFrequencies())
    throw Exception("NegBinomialFrequencySet::setFrequencies: wrong number of frequencies.");

  setFrequencies_(frequencies);
}

void NegBinomialFrequencySet::updateFrequencies_()
{
  const size_t n   = getNumberOfFrequencies();
  const double r   = getR();
  const double p   = r / (r + mu_);
  const double lp  = std::log(p);
  const double l1mp = std::log(1.0 - p);

  std::vector<double> freqs(n);

  for (size_t k = 0; k < n; ++k)
  {
    double logProb = std::lgamma(static_cast<double>(k) + r)
                   - std::lgamma(static_cast<double>(k) + 1.0)
                   - std::lgamma(r)
                   + r * lp
                   + static_cast<double>(k) * l1mp;
    freqs[k] = std::exp(logProb);
  }

  // Normalize over truncated support
  double sum = 0.0;
  for (double v : freqs)
    sum += v;

  if (sum <= 0.0)
    throw Exception("NegBinomialFrequencySet: normalization failed.");

  for (double& v : freqs)
    v /= sum;

  setFrequencies_(freqs);
}
