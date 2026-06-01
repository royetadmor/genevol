#include <iostream>
#include <string> 
#include <set>

// From bpp-core
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/Functions.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/App/BppApplication.h>

// From bpp-phyl
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>

// From bpp-seq
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// Local modules
#include "ModelParameters.h"
#include "LikelihoodUtils.h"
#include "MixtureModelLikelihoodFunction.h"
#include "GeneCountManager.h"
#include "TreeUtils.h"
#include "PoissonDistribution.h"
#include "WGDManager.h"

using namespace bpp;
using namespace std;

int main(int args, char **argv) {
    // Set model data and parameters
    BppApplication GenEvol(args, argv, "GenEvol");
    ModelParameters* m = new ModelParameters(GenEvol);

    // Get tree and rescale it
    Newick reader;
    std::shared_ptr<bpp::PhyloTree> tree_ = std::move(reader.readPhyloTree(m->treeFilePath_));
    double scale_tree_factor = TreeUtils::getTreeScalingFactor(m, tree_);
    tree_->scaleTree(scale_tree_factor);

    // Define substitution parameters
    auto paramMap = m->paramMap_;
    auto rateChangeType = m->rateChangeType_;
    auto constraintedParams = m->constraintedParams_;
    auto rDist = m->rDist_;

    // Calculate new likelihood
    auto likProc = LikelihoodUtils::createLikelihoodProcess(m, tree_, paramMap, rateChangeType, constraintedParams, rDist);
    std::cout << "Likelihood: " << likProc->getValue() << std::endl;
    if(std::isinf(likProc->getValue())) {
        std::cout << "Likelihood is inf, exiting" << std::endl;
        return 1;
    }
    // Optimization and assessment
    std::cout << "Starting optimization for new model" << std::endl;
    LikelihoodUtils::optimizeModelParametersOneDimension(likProc, m, m->optTolerance_, m->optNumIterations_);

    // Create mixture model, calculate likelihood and optimize
    if (m->useMixtureModel_) {
        auto geneCountManager = std::make_shared<GeneCountManager>(m, tree_);
        std::cout << "MM Likelihood: " << geneCountManager->getLikelihood() << std::endl;
        geneCountManager->optimizeMixtureModelParametersOneDimension(1e-4, 2);
        geneCountManager->printResults();
    }

    LikelihoodUtils::printResults(likProc, m->showRate4Site_);

    // WGD detection
    if (m->wgdThreshold_ > 0.0) {
        TreeUtils::printTopology(tree_);
        std::cout << "\nStarting WGD detection (threshold=" << m->wgdThreshold_ << ")" << std::endl;
        bpp::WGDManager wgdManager(m, tree_, likProc, m->wgdThreshold_);
        wgdManager.forwardPass();
        wgdManager.printResults();
        wgdManager.writeTree();
    }

    GenEvol.done();
    return 0;
}
