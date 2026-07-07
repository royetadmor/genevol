#ifndef GENEVOL_MODEL_ADEQUACY_UTILS_H
#define GENEVOL_MODEL_ADEQUACY_UTILS_H

#include <cmath>

#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

using namespace std;
namespace bpp {

class ModelAdequacyUtils {
public:
    ModelAdequacyUtils() {}
    virtual ~ModelAdequacyUtils() {}

    static double calculateAIC(SingleProcessPhyloLikelihood* lik, size_t extraParams = 0);

    // Chi-squared survival function: P(X > stat | df).
    static double chi2pvalue(double stat, int df);

private:
    // Series expansion of the regularized lower incomplete gamma P(a, x).
    static double gammaIncSeries(double a, double x);
};

} // namespace bpp

#endif // GENEVOL_MODEL_ADEQUACY_UTILS_H
