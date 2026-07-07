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


void WGDManager::optimizeParam(FunctionInterface* func, const std::string& paramName,
                               double lo, double hi, double tolerance)
{
    ParameterList params = func->getParameters();
    std::string actualName;
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName().find(paramName) != std::string::npos) {
            actualName = params[i].getName();
            break;
        }
    }
    if (actualName.empty())
        throw std::runtime_error("optimizeParam: parameter '" + paramName + "' not found.");

    auto f = std::shared_ptr<FunctionInterface>(func, [](FunctionInterface*) {});
    ExtendedBrentOptimizer optimizer(f);
    optimizer.setVerbose(0);
    optimizer.setProfiler(0);
    optimizer.setMessageHandler(0);
    optimizer.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer.setMaximumNumberOfEvaluations(100);
    optimizer.setBracketing(ExtendedBrentOptimizer::BRACKET_SIMPLE);
    optimizer.getStopCondition()->setTolerance(tolerance);
    optimizer.setInitialInterval(lo, hi);
    optimizer.init(params.createSubList(actualName));
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
    const std::string qParamName = "WGD_" + std::to_string(ins.wgdEdgeIdx) + ".q";
    auto tConstraint = dynamic_pointer_cast<IntervalConstraint>(
        posFunc.getParameters().getParameter("t")->getConstraint());
    for (int round = 0; round < 2; ++round) {
        optimizeParam(altLik, qParamName, 0.0, 1.0, 1e-4);
        optimizeParam(&posFunc, "t", tConstraint->getLowerBound(), tConstraint->getUpperBound(), m_->optTolerance_);
    }

    // t adds one extra parameter not counted by calculateAIC, so penalize by +2
    double deltaAIC = baseAIC - (ModelAdequacyUtils::calculateAIC(altLik, acceptedTCount_) + 2.0);
    double q = 0.5;
    ParameterList ps = altLik->getParameters();
    for (size_t pi = 0; pi < ps.size(); ++pi) {
        if (ps[pi].getName().find(qParamName) != std::string::npos) {
            q = ps[pi].getValue(); break;
        }
    }

    best.childId  = childId;
    best.deltaAIC = deltaAIC;
    best.q        = q;
    best.t        = posFunc.getParameterValue("t");
    best.lik      = altLik;

    TreeUtils::removeWGDNode(tree_, ins, nextEdgeIdx_);

    std::cout << "  Branch to node " << childId
              << ": best t=" << best.t << "  ΔAIC=" << best.deltaAIC << "  q=" << best.q << std::endl;

    return best;
}

uint WGDManager::getLeafCount(uint nodeId) const
{
    auto node = tree_->getNode(nodeId);
    if (tree_->isLeaf(node))
        return 1;
    uint count = 0;
    for (const auto& child : tree_->getSons(node))
        count += getLeafCount(tree_->getNodeIndex(child));
    return count;
}

void WGDManager::forwardPass()
{
    double baseAIC = ModelAdequacyUtils::calculateAIC(baseLik_);
    std::cout << "WGD forward pass — baseline AIC: " << baseAIC << std::endl;

    while (true) {
        auto currentParams = extractRateParams(baseLik_);
        auto currentRDist  = extractRDist(baseLik_);

        // Track the best candidate incrementally: shallowest depth, then highest ΔAIC
        CandidateResult best;

        for (uint childId : getCandidates()) {
            CandidateResult res = evaluateCandidate(childId, baseAIC, currentParams, currentRDist);
            res.leafCount = getLeafCount(childId);

            if (res.deltaAIC <= threshold_) {
                if (res.lik) LikelihoodUtils::deleteLikelihoodProcess(res.lik);
                continue;
            }

            if (best.lik == nullptr || beats(res, best)) {
                if (best.lik) LikelihoodUtils::deleteLikelihoodProcess(best.lik);
                best = std::move(res);
            } else {
                LikelihoodUtils::deleteLikelihoodProcess(res.lik);
            }
        }

        if (best.lik == nullptr) {
            std::cout << "WGD forward pass complete. Found " << results_.size() << " duplications." << std::endl;
            break;
        }

        auto bestChild = tree_->getNode(best.childId);
        auto acceptedIns = TreeUtils::insertWGDNode(tree_, bestChild, nextNodeIdx_, nextEdgeIdx_, best.t);

        std::cout << "Accepted WGD on branch to node " << best.childId
                  << "  t=" << best.t << "  ΔAIC=" << best.deltaAIC << "  q=" << best.q << std::endl;

        wgdQMap_[acceptedIns.wgdEdgeIdx] = best.q;

        ownedLiks_.push_back(best.lik);
        baseLik_ = best.lik;
        acceptedTCount_++;
        baseAIC = ModelAdequacyUtils::calculateAIC(baseLik_, acceptedTCount_);

        WGDResult res;
        res.childNodeId = best.childId;
        res.wgdEdgeIdx  = acceptedIns.wgdEdgeIdx;
        res.q           = best.q;
        res.t           = best.t;
        res.deltaAIC    = best.deltaAIC;
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

void WGDManager::testWGD()
{
    // Collect zero-length edge IDs (user-specified WGD events)
    std::vector<uint> wgdEdgeIds;
    for (auto& node : tree_->getAllNodes()) {
        if (tree_->getNodeIndex(node) == tree_->getRootIndex()) continue;
        auto edge = tree_->getEdgeToFather(node);
        if (edge->getLength() == 0.0)
            wgdEdgeIds.push_back(tree_->getEdgeIndex(edge));
    }

    if (wgdEdgeIds.empty()) {
        return;
    }

    int k = static_cast<int>(wgdEdgeIds.size());
    std::cout << "\nWGD test mode: " << k << " event(s) found in input tree." << std::endl;

    // Build null model: all q = 0, re-optimize rate params only
    std::map<uint, double> nullQMap;
    for (uint edgeId : wgdEdgeIds)
        nullQMap[edgeId] = 0.0;

    auto nullLik = LikelihoodUtils::createLikelihoodProcess(
        m_, tree_, m_->paramMap_, m_->rateChangeType_, m_->constraintedParams_, m_->rDist_, nullQMap, 0.0);

    m_->fixedParams_.push_back("WGD");
    LikelihoodUtils::optimizeModelParametersOneDimension(nullLik, m_, m_->optTolerance_, m_->optNumIterations_);
    m_->fixedParams_.pop_back();

    // Print per-event q values from alt model
    std::cout << "\n=== WGD Test Results ===" << std::endl;
    ParameterList altParams = baseLik_->getParameters();
    for (uint edgeId : wgdEdgeIds) {
        std::string qName = "WGD_" + std::to_string(edgeId) + ".q";
        double q = -1.0;
        for (size_t i = 0; i < altParams.size(); ++i) {
            if (altParams[i].getName().find(qName) != std::string::npos) {
                q = altParams[i].getValue(); break;
            }
        }
        std::cout << "  Edge " << edgeId << ": q = " << q << std::endl;
    }

    if (m_->modelCriterion_ == "AIC") {
        double altAIC  = ModelAdequacyUtils::calculateAIC(baseLik_);
        double nullAIC = ModelAdequacyUtils::calculateAIC(nullLik);
        double deltaAIC = nullAIC - altAIC;
        std::cout << "  Alt  AIC=" << altAIC  << std::endl;
        std::cout << "  Null AIC=" << nullAIC << std::endl;
        std::cout << "  ΔAIC (null - alt) = " << deltaAIC << std::endl;
        std::cout << "  Decision (AIC, threshold=" << m_->wgdThreshold_ << "): "
                  << (deltaAIC > m_->wgdThreshold_ ? "WGD supported" : "WGD not supported") << std::endl;
    } else {
        double lrt  = 2.0 * (nullLik->getValue() - baseLik_->getValue());
        double pval = ModelAdequacyUtils::chi2pvalue(lrt, k);
        std::cout << "  LRT = " << lrt << "  df = " << k << std::endl;
        std::cout << "  p-value = " << pval << std::endl;
        std::cout << "  Decision (LRT, α=0.05): " << (pval < 0.05 ? "WGD supported" : "WGD not supported") << std::endl;
    }

    LikelihoodUtils::deleteLikelihoodProcess(nullLik);
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

