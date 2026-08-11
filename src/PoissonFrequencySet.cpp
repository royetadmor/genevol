// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PoissonFrequencySet.h"

#include <Bpp/Numeric/NumTools.h>
#include <cmath>
#include <stdexcept>

using namespace bpp;
using namespace std;

PoissonFrequencySet::PoissonFrequencySet(
  shared_ptr<const StateMapInterface> stateMap,
  double lambda,
  bool optimize)
  :
  AbstractFrequencySet(stateMap, "Poisson.", "Poisson"),
  lambda_(lambda),
  optimize_(optimize)
{
  if (lambda_ <= 0.0)
    throw Exception("PoissonFrequencySet: lambda must be > 0.");

  if (optimize_)
    addParameter_(new Parameter("Poisson.lambda", lambda_, Parameter::R_PLUS_STAR));

  updateFrequencies_();
}

void PoissonFrequencySet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequencySet::fireParameterChanged(parameters);
  if (optimize_)
    lambda_ = getParameterValue("lambda");
  updateFrequencies_();
}

void PoissonFrequencySet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getNumberOfFrequencies())
    throw Exception("PoissonFrequencySet::setFrequencies: wrong number of frequencies.");

  setFrequencies_(frequencies);
}

void PoissonFrequencySet::updateFrequencies_()
{
  const size_t n = getNumberOfFrequencies();
  vector<double> freqs(n);

  double logLambda = std::log(lambda_);
  double logFactorial = 0.0;

  for (size_t k = 0; k < n; ++k)
  {
    if (k > 0)
      logFactorial += std::log(static_cast<double>(k));
    freqs[k] = std::exp(-lambda_ + static_cast<double>(k) * logLambda - logFactorial);
  }

  double sum = 0.0;
  for (double v : freqs)
    sum += v;

  if (sum <= 0.0)
    throw Exception("PoissonFrequencySet: normalization failed.");

  for (double& v : freqs)
    v /= sum;

  setFrequencies_(freqs);
}
