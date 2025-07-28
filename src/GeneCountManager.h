#ifndef GENEVOL_GENECOUNTMANAGER_H
#define GENEVOL_GENECOUNTMANAGER_H

#include <memory>

#include <Bpp/Phyl/Tree/PhyloTree.h>


#include "ModelParameters.h"
#include "MixtureModelLikelihoodFunction.h" // for MixtureModelLikelihoodFunction
#include "LikelihoodUtils.h"

namespace bpp {

class GeneCountManager {
public:
    // Constructor
    GeneCountManager(ModelParameters* m, PhyloTree* tree) {
        m_ = m;
        tree_ = tree;
        likelihoodFunction_ = std::make_shared<MixtureModelLikelihoodFunction>(m, tree_);
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
    std::shared_ptr<MixtureModelLikelihoodFunction> likelihoodFunction_;
};

} // namespace bpp

#endif
