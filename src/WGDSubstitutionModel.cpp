#include "WGDSubstitutionModel.h"

using namespace bpp;
using namespace std;

WGDSubstitutionModel::WGDSubstitutionModel(
    std::shared_ptr<GeneCountSubstitutionModel> baseModel,
    double q,
    int maxState)
    : AbstractParameterAliasable("WGD."),
      AbstractSubstitutionModel(baseModel->getAlphabet(), baseModel->getStateMap(), "WGD."),
      W_(baseModel->getNumberOfStates(), baseModel->getNumberOfStates()),
      maxState_(maxState)
{
    addParameter_(new Parameter("WGD.q", q, make_shared<IntervalConstraint>(0, 1.0, true, true)));
    isScalable_ = false;
    computeWMatrix_();
}

void WGDSubstitutionModel::computeWMatrix_() {
    const double q = getQ();
    const double logQ = log(q);
    const double log1mQ = log(1.0 - q);
    const size_t n = size_;

    MatrixTools::fill(W_, 0.0);

    for (size_t i = 0; i < n; i++) {
        if (i == 0) {
            W_(0, 0) = 1.0;
            continue;
        }

        double rowSum = 0.0;
        size_t jMax = std::min(2 * i, n - 1);  // cap at maxState

        for (size_t j = i; j <= jMax; j++) {
            size_t retained = j - i;          // number of duplicate copies kept
            size_t lost     = i - retained;   // = 2i - j

            // log C(i, j-i) + (j-i)*log(q) + (2i-j)*log(1-q)
            double logProb = lgamma(i + 1) - lgamma(retained + 1) - lgamma(lost + 1);
            if (retained > 0) logProb += retained * logQ;
            if (lost     > 0) logProb += lost     * log1mQ;

            W_(i, j) = exp(logProb);
            if (std::isnan(W_(i, j)) || std::isinf(W_(i, j)))
                W_(i, j) = 0.0;
            rowSum += W_(i, j);
        }

        // Renormalize if truncation at maxState cut off part of the distribution
        if (rowSum > 0.0 && rowSum < 1.0 - 1e-10) {
            for (size_t j = i; j <= jMax; j++)
                W_(i, j) /= rowSum;
        }
    }
}

void WGDSubstitutionModel::updateMatrices_() {
    computeWMatrix_();
}

const Matrix<double>& WGDSubstitutionModel::getPij_t(double t) const {
    // WGD is an instantaneous event — branch length is irrelevant.
    // Simply return the WGD transition matrix W(q).
    MatrixTools::copy(W_, pijt_);

    return pijt_;
}

const Matrix<double>& WGDSubstitutionModel::getdPij_dt(double t) const {
    // WGD is instantaneous — derivative with respect to branch length is zero.
    MatrixTools::fill(dpijt_, 0.0);
    return dpijt_;
}

const Matrix<double>& WGDSubstitutionModel::getd2Pij_dt2(double t) const {
    // WGD is instantaneous — second derivative with respect to branch length is zero.
    MatrixTools::fill(d2pijt_, 0.0);
    return d2pijt_;
}
