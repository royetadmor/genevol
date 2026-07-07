#include "ModelAdequacyUtils.h"

using namespace bpp;
using namespace std;

double ModelAdequacyUtils::calculateAIC(SingleProcessPhyloLikelihood* lik, size_t extraParams)
{
    ParameterList params = lik->getSubstitutionProcess()->getIndependentParameters();
    size_t numOfParams = 0;
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName().find("BrLen") == std::string::npos)
            numOfParams++;
    }
    numOfParams += extraParams;
    return 2.0 * lik->getValue() + 2.0 * numOfParams;
}

double ModelAdequacyUtils::gammaIncSeries(double a, double x)
{
    double sum = 1.0 / a, term = 1.0 / a;
    for (int n = 1; n <= 500; n++) {
        term *= x / (a + n);
        sum += term;
        if (std::abs(term) < 1e-14 * std::abs(sum)) break;
    }
    return std::exp(-x + a * std::log(x) - std::lgamma(a)) * sum;
}

double ModelAdequacyUtils::chi2pvalue(double stat, int df)
{
    if (stat <= 0.0) return 1.0;
    return 1.0 - gammaIncSeries(df / 2.0, stat / 2.0);
}

