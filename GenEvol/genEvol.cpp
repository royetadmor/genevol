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
    m->validateTree(tree_);
    double scale_tree_factor = TreeUtils::getTreeScalingFactor(m, tree_);
    tree_->scaleTree(scale_tree_factor);

    // Define substitution parameters
    auto paramMap = m->paramMap_;
    auto rateChangeType = m->rateChangeType_;
    auto constraintedParams = m->constraintedParams_;
    auto rDist = m->rDist_;

    // Multi-start optimization
    auto likProc = LikelihoodUtils::multiStartOptimize(m, tree_, rateChangeType, constraintedParams, rDist);
    if (!likProc) {
        std::cout << "All starting points gave infinite likelihood, exiting" << std::endl;
        return 1;
    }

    // Create mixture model, calculate likelihood and optimize
    if (m->useMixtureModel_) {
        auto geneCountManager = std::make_shared<GeneCountManager>(m, tree_);
        std::cout << "MM Likelihood: " << geneCountManager->getLikelihood() << std::endl;
        geneCountManager->optimizeMixtureModelParametersOneDimension(1e-4, 2);
        geneCountManager->printResults();
    }

    LikelihoodUtils::printResults(likProc, m->showRate4Site_);

    // WGD analysis
    if (m->wgdMode_ == "detect") {
        TreeUtils::printTopology(tree_);
        std::cout << "\nStarting WGD detection (threshold=" << m->wgdThreshold_ << ")" << std::endl;
        bpp::WGDManager wgdManager(m, tree_, likProc, m->wgdThreshold_);
        wgdManager.forwardPass();
        wgdManager.printDetectionResults();
        TreeUtils::writeTree(tree_, "wgd_tree.nwk");
    } else if (m->wgdMode_ == "test") {
        std::vector<uint> wgdEdgeIds = TreeUtils::collectWgdEdgesInOrder(tree_);
        std::map<uint, double> fixedEdges;
        std::vector<uint> freeEdgeIds;
        for (size_t i = 0; i < wgdEdgeIds.size(); ++i) {
            if (i < m->fixedWgdQ_.size())
                fixedEdges[wgdEdgeIds[i]] = m->fixedWgdQ_[i];
            else
                freeEdgeIds.push_back(wgdEdgeIds[i]);
        }

        SingleProcessPhyloLikelihood* baseLik = likProc;
        // Assiging all fixed WGDs (if any) and reoptimizing parameters
        if (!fixedEdges.empty()) {
            baseLik = LikelihoodUtils::createLikelihoodProcess(
                m, tree_, paramMap, rateChangeType, constraintedParams, rDist, fixedEdges, m->rootLambda_);
            for (const auto& kv : fixedEdges) {
                m->fixedParams_.push_back("WGD_" + std::to_string(kv.first) + ".q");
            }
            LikelihoodUtils::optimizeModelParametersOneDimension(baseLik, m, m->optTolerance_, m->optNumIterations_);
            for (size_t i = 0; i < fixedEdges.size(); ++i) {
                m->fixedParams_.pop_back();
            }
        }

        // If there are free q values, we optimize them here.
        // Otherwise, we just print the results and exit
        WGDManager wgdManager(m, tree_, baseLik, m->wgdThreshold_);
        wgdManager.testWGD(fixedEdges, freeEdgeIds);
    }

    GenEvol.done();
    return 0;
}
