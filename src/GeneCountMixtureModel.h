#ifndef GENEVOL_GENECOUNTMIXTUREMODEL_H
#define GENEVOL_GENECOUNTMIXTUREMODEL_H

#include <memory>

#include <Bpp/Phyl/Tree/PhyloTree.h>


#include "ModelParameters.h"
#include "MixtureModelLikelihoodFunction.h" // for AlphaLikelihoodFunction
#include "LikelihoodUtils.h"

namespace bpp {

class GeneCountMixtureModel {
public:
    // Constructor
    GeneCountMixtureModel(ModelParameters* m, PhyloTree* tree) {
        m_ = m;
        tree_ = tree;
        likelihoodFunction_ = std::make_shared<AlphaLikelihoodFunction>(m, tree_);
        vectorOfLikelihoods = likelihoodFunction_->getLikelihoodProcesses();
    }

    // Optimize mixture model parameters
    void optimizeMixtureModelParametersOneDimension(double tol, unsigned int maxNumOfIterations, bool mixed=false, unsigned curentIterNum=0);

    // Get final log likelihood
    double getLikelihood() const { return likelihoodFunction_->getValue(); }
    double getParameterValueByName(string name) const { return likelihoodFunction_->getParameterValueByName(name); }

private:
    ModelParameters* m_;
    PhyloTree* tree_;
    std::vector<SingleProcessPhyloLikelihood*> vectorOfLikelihoods;
    std::shared_ptr<AlphaLikelihoodFunction> likelihoodFunction_;
};

} // namespace bpp

#endif
