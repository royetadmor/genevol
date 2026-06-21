#ifndef GENEVOL_WGD_POSITION_FUNCTION_H
#define GENEVOL_WGD_POSITION_FUNCTION_H

#include <iostream>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

namespace bpp {

/**
 * Wraps a SingleProcessPhyloLikelihood so that Brent can optimize the WGD
 * position `t ∈ (0,1)` along a branch.  When `t` changes, the upper and lower
 * branch lengths are updated via BrLen{id} parameters on the likelihood before
 * re-evaluation.
 *
 * The likelihood is not owned by this class.
 */
class WGDPositionFunction :
    public FunctionInterface,
    public AbstractParameterAliasable
{
private:
    SingleProcessPhyloLikelihood* lik_;
    std::string upperBrLenName_;
    std::string lowerBrLenName_;
    double origLen_;
    double fval_;

public:
    WGDPositionFunction(SingleProcessPhyloLikelihood* lik,
                        uint upperBranchId,
                        uint lowerBranchId,
                        double origLen,
                        double initT)
        : AbstractParameterAliasable(""),
          FunctionInterface(),
          lik_(lik),
          upperBrLenName_("BrLen" + std::to_string(upperBranchId)),
          lowerBrLenName_("BrLen" + std::to_string(lowerBranchId)),
          origLen_(origLen),
          fval_(lik->getValue())
    {
        addParameter_(new Parameter("t", initT,
            std::make_shared<IntervalConstraint>(1e-6, 1.0 - 1e-6, false, false)));
        fireParameterChanged(getParameters());
    }

    Clonable* clone() const override { return new WGDPositionFunction(*this); }

    double getValue() const override { return fval_; }

    void setParameters(const ParameterList& parameters) override {
        matchParametersValues(parameters);
    }

    void fireParameterChanged(const ParameterList&) override {
        double t = getParameterValue("t");
        lik_->getParameters().getParameter(upperBrLenName_)->setValue(origLen_ * t);
        lik_->getParameters().getParameter(lowerBrLenName_)->setValue(origLen_ * (1.0 - t));
        fval_ = lik_->getValue();
    }
};

} // namespace bpp

#endif // GENEVOL_WGD_POSITION_FUNCTION_H
