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


void WGDManager::optimizeQ(SingleProcessPhyloLikelihood* lik, uint newEdgeIdx)
{
    const std::string targetParam = "WGD_" + std::to_string(newEdgeIdx) + ".q";

    ParameterList params = lik->getParameters();
    bool found = false;
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName().find(targetParam) != std::string::npos) {
            const std::string actualParam = params[i].getName();
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
            optimizer.init(params.createSubList(actualParam));
            optimizer.optimize();
            found = true;
            break;
        }
    }
    if (!found)
        throw std::runtime_error("optimizeQ: parameter '" + targetParam + "' not found in likelihood.");
}

void WGDManager::optimizeT(WGDPositionFunction* posFunc)
{
    ParameterList params = posFunc->getParameters();
    std::string tParamName;
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName() == "t") {
            tParamName = params[i].getName();
            break;
        }
    }
    if (tParamName.empty())
        throw std::runtime_error("optimizeT: 't' parameter not found in WGDPositionFunction.");

    auto constraint = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(tParamName)->getConstraint());
    double tMin = constraint->getLowerBound();
    double tMax = constraint->getUpperBound();
    auto f = std::shared_ptr<FunctionInterface>(posFunc, [](FunctionInterface*) {});
    ExtendedBrentOptimizer optimizer(f);
    optimizer.setVerbose(0);
    optimizer.setProfiler(0);
    optimizer.setMessageHandler(0);
    optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer.setMaximumNumberOfEvaluations(100);
    optimizer.setBracketing(ExtendedBrentOptimizer::BRACKET_SIMPLE);
    optimizer.getStopCondition()->setTolerance(m_->optTolerance_);
    optimizer.setInitialInterval(tMin, tMax);
    optimizer.init(params.createSubList(tParamName));
    optimizer.optimize();
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

WGDManager::CandidateResult WGDManager::evaluateCandidate(
    uint childId, double baseAIC,
    const std::map<int, std::vector<double>>& currentParams,
    std::shared_ptr<DiscreteDistributionInterface> currentRDist)
{
    auto candChild = tree_->getNode(childId);
    CandidateResult best;

    auto ins = TreeUtils::insertWGDNode(tree_, candChild, nextNodeIdx_, nextEdgeIdx_, 0.5);

    auto altLik = LikelihoodUtils::createLikelihoodProcess(
        m_, tree_, currentParams, m_->rateChangeType_,
        m_->constraintedParams_, currentRDist, wgdQMap_, 0.5);

    uint upperBranchId = tree_->getEdgeIndex(tree_->getEdgeToFather(ins.wgdUpper));
    uint lowerBranchId = tree_->getEdgeIndex(tree_->getEdgeToFather(candChild));

    // Alternating optimization: 2 rounds of (q, t)
    WGDPositionFunction posFunc(altLik, upperBranchId, lowerBranchId, ins.origLen, 0.5);
    for (int round = 0; round < 2; ++round) {
        optimizeQ(altLik, ins.wgdEdgeIdx);
        optimizeT(&posFunc);
    }

    // t adds one extra parameter not counted by calculateAIC, so penalize by +2
    double deltaAIC = baseAIC - (LikelihoodUtils::calculateAIC(altLik, acceptedTCount_) + 2.0);
    const std::string qParamName = "WGD_" + std::to_string(ins.wgdEdgeIdx) + ".q";
    double q = 0.5;
    ParameterList ps = altLik->getParameters();
    for (size_t pi = 0; pi < ps.size(); ++pi) {
        if (ps[pi].getName().find(qParamName) != std::string::npos) {
            q = ps[pi].getValue(); break;
        }
    }

    best.deltaAIC = deltaAIC;
    best.q        = q;
    best.t        = posFunc.getParameterValue("t");
    best.lik      = altLik;

    TreeUtils::removeWGDNode(tree_, ins, nextEdgeIdx_);

    std::cout << "  Branch to node " << childId
              << ": best t=" << best.t << "  ΔAIC=" << best.deltaAIC << "  q=" << best.q << std::endl;

    return best;
}

void WGDManager::forwardPass()
{
    double baseAIC = LikelihoodUtils::calculateAIC(baseLik_);
    std::cout << "WGD forward pass — baseline AIC: " << baseAIC << std::endl;

    while (true) {
        double bestDeltaAIC = 0.0;
        int    bestChildId  = -1;
        double bestQ        = 0.5;
        double bestT        = 0.5;
        SingleProcessPhyloLikelihood* bestLik = nullptr;

        auto currentParams = extractRateParams(baseLik_);
        auto currentRDist  = extractRDist(baseLik_);

        for (uint childId : getCandidates()) {
            CandidateResult res = evaluateCandidate(childId, baseAIC, currentParams, currentRDist);

            if (res.deltaAIC > bestDeltaAIC) {
                if (bestLik) LikelihoodUtils::deleteLikelihoodProcess(bestLik);
                bestDeltaAIC = res.deltaAIC;
                bestChildId  = childId;
                bestQ        = res.q;
                bestT        = res.t;
                bestLik      = res.lik;
            } else {
                if (res.lik) LikelihoodUtils::deleteLikelihoodProcess(res.lik);
            }
        }

        if (bestDeltaAIC <= threshold_ || bestChildId == -1) {
            if (bestLik) LikelihoodUtils::deleteLikelihoodProcess(bestLik);
            std::cout << "WGD forward pass complete. Found " << results_.size() << " duplications." << std::endl;
            break;
        }

        auto bestChild = tree_->getNode(static_cast<uint>(bestChildId));
        auto acceptedIns = TreeUtils::insertWGDNode(tree_, bestChild, nextNodeIdx_, nextEdgeIdx_, bestT);

        std::cout << "Accepted WGD on branch to node " << bestChildId
                  << "  t=" << bestT << "  ΔAIC=" << bestDeltaAIC << "  q=" << bestQ << std::endl;

        wgdQMap_[acceptedIns.wgdEdgeIdx] = bestQ;

        ownedLiks_.push_back(bestLik);
        baseLik_ = bestLik;
        acceptedTCount_++;
        baseAIC = LikelihoodUtils::calculateAIC(baseLik_, acceptedTCount_);

        WGDResult res;
        res.childNodeId  = static_cast<uint>(bestChildId);
        res.wgdEdgeIdx   = acceptedIns.wgdEdgeIdx;
        res.q            = bestQ;
        res.t            = bestT;
        res.deltaAIC     = bestDeltaAIC;
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
    std::cout << "  #   Child-node       q        ΔAIC" << std::endl;
    for (size_t i = 0; i < results_.size(); ++i) {
        const auto& r = results_[i];
        std::cout << "  " << (i + 1)
                  << "       " << r.childNodeId
                  << "      " << r.q
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

