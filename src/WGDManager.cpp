#include "WGDManager.h"
#include "ExtendedBrentOptimizer.h"
#include "TreeUtils.h"

#include <Bpp/Numeric/AutoParameter.h>

using namespace bpp;
using namespace std;

WGDManager::WGDManager(ModelParameters* m,
                       std::shared_ptr<PhyloTree> tree,
                       SingleProcessPhyloLikelihood* baseLik,
                       double threshold)
    : m_(m), tree_(tree), baseLik_(baseLik), threshold_(threshold)
{
    // Seed nextEdgeIdx_ and nextNodeIdx_ from the tree's current max so we never collide
    for (auto& n : tree->getAllNodes()) {
        if (!tree->hasNodeIndex(n)) continue;
        nextNodeIdx_ = std::max(nextNodeIdx_, tree->getNodeIndex(n) + 1);
        if (tree->getNodeIndex(n) == tree->getRootIndex()) continue;
        auto e = tree->getEdgeToFather(n);
        if (tree->hasEdgeIndex(e))
            nextEdgeIdx_ = std::max(nextEdgeIdx_, tree->getEdgeIndex(e) + 1);
    }
}


void WGDManager::optimizeQ(SingleProcessPhyloLikelihood* lik)
{
    auto f = std::shared_ptr<SecondOrderDerivable>(lik, [](SecondOrderDerivable*) {});
    ExtendedBrentOptimizer optimizer(f);
    optimizer.setVerbose(0);
    optimizer.setProfiler(0);
    optimizer.setMessageHandler(0);
    optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer.setMaximumNumberOfEvaluations(100);
    optimizer.setBracketing(ExtendedBrentOptimizer::BRACKET_SIMPLE);
    optimizer.getStopCondition()->setTolerance(1e-4);
    optimizer.setInitialInterval(0, 1.0);

    ParameterList params = lik->getParameters();
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName().find("WGD.q") != std::string::npos) {
            optimizer.init(params.createSubList(params[i].getName()));
            optimizer.optimize();
            break;
        }
    }
}

std::map<int, std::vector<double>> WGDManager::extractRateParams(SingleProcessPhyloLikelihood* lik) const
{
    static const std::vector<std::pair<std::string, int>> nameToType = {
        {"loss",        GeneCountSubstitutionModel::paramType::LOSS},
        {"gain",        GeneCountSubstitutionModel::paramType::GAIN},
        {"innovation",  GeneCountSubstitutionModel::paramType::INNOVATION},
        {"elimination", GeneCountSubstitutionModel::paramType::ELIMINATION}
    };

    std::map<int, std::vector<double>> result = m_->paramMap_;

    ParameterList params = lik->getParameters();
    for (size_t i = 0; i < params.size(); ++i) {
        const std::string& name = params[i].getName();
        for (const auto& kv : nameToType) {
            if (name.find(kv.first) != std::string::npos) {
                int idx = LikelihoodUtils::getParamIndex(name);
                auto& vec = result[kv.second];
                if (idx >= static_cast<int>(vec.size()))
                    vec.resize(idx + 1, 0.0);
                vec[idx] = params[i].getValue();
                break;
            }
        }
    }
    return result;
}

std::shared_ptr<DiscreteDistributionInterface> WGDManager::extractRDist(SingleProcessPhyloLikelihood* lik) const
{
    size_t nCat = m_->rDist_->getNumberOfCategories();
    if (nCat == 1) {
        return m_->rDist_;
    }

    double alpha = 1.0;
    ParameterList likParams = lik->getParameters();
    for (size_t j = 0; j < likParams.size(); ++j) {
        if (likParams[j].getName().find("alpha") != std::string::npos) {
            alpha = likParams[j].getValue();
            break;
        }
    }
    return std::make_shared<GammaDiscreteRateDistribution>(nCat, alpha);
}

void WGDManager::forwardPass()
{
    double baseAIC = LikelihoodUtils::calculateAIC(baseLik_);
    std::cout << "WGD forward pass — baseline AIC: " << baseAIC << std::endl;

    while (true) {
        double bestLRT      = 0.0;
        uint   bestChildId  = 0;
        double bestQ        = 0.5;
        SingleProcessPhyloLikelihood* bestLik = nullptr;

        auto currentParams = extractRateParams(baseLik_);
        auto currentRDist  = extractRDist(baseLik_);
        auto candidateIds  = getCandidates();

        for (uint childId : candidateIds) {
            auto candChild = tree_->getNode(childId);
            auto [wgdNode, origLen] = TreeUtils::insertWGDNode(tree_, candChild, nextNodeIdx_, nextEdgeIdx_);

            // Null: q fixed at 0
            auto nullLik = LikelihoodUtils::createLikelihoodProcess(
                m_, tree_, currentParams, m_->rateChangeType_,
                m_->constraintedParams_, currentRDist, 0.0);
            
            // Alternative: q optimized
            auto altLik = LikelihoodUtils::createLikelihoodProcess(
                m_, tree_, currentParams, m_->rateChangeType_,
                m_->constraintedParams_, currentRDist, 0.5);
            optimizeQ(altLik);

            // LRT = 2 * (logL_alt - logL_null); getValue() returns -logL
            double lrt      = 2.0 * (nullLik->getValue() - altLik->getValue());
            double deltaAIC = baseAIC - LikelihoodUtils::calculateAIC(altLik);

            LikelihoodUtils::deleteLikelihoodProcess(nullLik);

            double q = 0.5;
            ParameterList ps = altLik->getParameters();
            for (size_t pi = 0; pi < ps.size(); ++pi) {
                if (ps[pi].getName().find("WGD.q") != std::string::npos) {
                    q = ps[pi].getValue(); break;
                }
            }

            std::cout << "  Branch to node " << childId
                      << ": LRT=" << lrt << " ΔAIC=" << deltaAIC << " q=" << q << std::endl;

            if (lrt > bestLRT) {
                if (bestLik) LikelihoodUtils::deleteLikelihoodProcess(bestLik);
                bestLRT     = lrt;
                bestChildId = childId;
                bestQ       = q;
                bestLik     = altLik;
            } else {
                LikelihoodUtils::deleteLikelihoodProcess(altLik);
            }

            TreeUtils::removeWGDNode(tree_, wgdNode, origLen, nextEdgeIdx_);
        }

        if (bestLRT <= threshold_ || bestChildId == 0) {
            if (bestLik) LikelihoodUtils::deleteLikelihoodProcess(bestLik);
            std::cout << "WGD forward pass complete. Found " << results_.size() << " duplications." << std::endl;
            break;
        }

        // Accept best WGD: insert permanently into tree_
        auto bestChild = tree_->getNode(bestChildId);
        TreeUtils::insertWGDNode(tree_, bestChild, nextNodeIdx_, nextEdgeIdx_);

        double bestDeltaAIC = baseAIC - LikelihoodUtils::calculateAIC(bestLik);
        std::cout << "Accepted WGD on branch to node " << bestChildId
                  << "  LRT=" << bestLRT << "  ΔAIC=" << bestDeltaAIC << "  q=" << bestQ << std::endl;

        baseAIC = LikelihoodUtils::calculateAIC(bestLik);
        ownedLiks_.push_back(bestLik);
        baseLik_ = bestLik;

        WGDResult res;
        res.childNodeId = bestChildId;
        res.q           = bestQ;
        res.lrt         = bestLRT;
        res.deltaAIC    = bestDeltaAIC;
        results_.push_back(res);
    }
}

std::vector<uint> WGDManager::getCandidates() const
{
    std::vector<uint> candidateIds;
    uint rootId = tree_->getRootIndex();
    for (auto& node : tree_->getAllNodes()) {
        uint nid = tree_->getNodeIndex(node);
        if (nid == rootId) continue;
        auto edge = tree_->getEdgeToFather(node);
        if (edge->getLength() <= 0.0) continue;
        auto parent = tree_->getFatherOfNode(node);
        if (parent->hasName() && parent->getName() == "WGD") continue;
        candidateIds.push_back(nid);
    }
    return candidateIds;
}

void WGDManager::printResults() const
{
    if (results_.empty()) {
        std::cout << "No WGD events detected." << std::endl;
        return;
    }
    std::cout << "\n=== WGD Detection Results ===" << std::endl;
    std::cout << "  #   Child-node       q        LRT      ΔAIC" << std::endl;
    for (size_t i = 0; i < results_.size(); ++i) {
        const auto& r = results_[i];
        std::cout << "  " << (i + 1)
                  << "       " << r.childNodeId
                  << "      " << r.q
                  << "   " << r.lrt
                  << "   " << r.deltaAIC << std::endl;
    }
    TreeUtils::printTopology(tree_);
}

void WGDManager::writeTree(const std::string& outputPath) const
{
    bpp::Newick newick;
    std::cout << "Newick tree: ";
    newick.writePhyloTree(*tree_, std::cout);
    std::cout << std::endl;

    std::ofstream out(outputPath);
    if (!out)
        throw std::runtime_error("Cannot open output file: " + outputPath);
    newick.writePhyloTree(*tree_, out);
    std::cout << "Tree with WGD events written to: " << outputPath << std::endl;
}

