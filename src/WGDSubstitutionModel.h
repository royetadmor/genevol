#ifndef GENEVOL_WGD_SUBSTITUTION_MODEL_H
#define GENEVOL_WGD_SUBSTITUTION_MODEL_H

#include <cmath>
#include <memory>

#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include "GeneCountSubstitutionModel.h"

namespace bpp {

/**
 * WGD substitution model based on Rabbier (2014).
 *  Tranistion probability is computed as follows:
 *  P(j|i) = C(i, j-i) * q^(j-i) * (1-q)^(2i-j),  for i <= j <= 2i
 *  where q = retention probability of each duplicate copy.
 *
 * W(j|k) is the WGD transition matrix: a family of size k doubles to 2k,
 * then each gene is independently lost with probability q (Rabbier 2014).
 * W(j|k) = Binomial(2k, j, 1-q), truncated and renormalized at maxState.
 *
 * The only free parameter is q in (0, 1).
 */
class WGDSubstitutionModel : public AbstractSubstitutionModel {
private:
    mutable RowMatrix<double> W_;
    int maxState_;

public:
    WGDSubstitutionModel(std::shared_ptr<GeneCountSubstitutionModel> baseModel, double q, int maxState);
    WGDSubstitutionModel(const WGDSubstitutionModel& model) = default;
    WGDSubstitutionModel& operator=(const WGDSubstitutionModel&) = default;
    virtual ~WGDSubstitutionModel() {}

    WGDSubstitutionModel* clone() const override { return new WGDSubstitutionModel(*this); }
    std::string getName() const override { return "WGD"; }

    const Matrix<double>& getPij_t    (double t) const override;
    const Matrix<double>& getdPij_dt  (double t) const override;
    const Matrix<double>& getd2Pij_dt2(double t) const override;

    double getQ() const { return getParameterValue("q"); }

protected:
    void updateMatrices_() override;

private:
    void computeWMatrix_();
};

} // namespace bpp

#endif // GENEVOL_WGD_SUBSTITUTION_MODEL_H
