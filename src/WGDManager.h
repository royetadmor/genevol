#ifndef GENEVOL_WGD_MANAGER_H
#define GENEVOL_WGD_MANAGER_H

#include <vector>
#include <map>
#include <memory>
#include <utility>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

#include "ModelParameters.h"
#include "GeneCountSubstitutionModel.h"
#include "LikelihoodUtils.h"

namespace bpp {

/**
 * WGDManager — greedy forward-pass WGD placement.
 *
 * Receives an already-optimized baseline likelihood.  For each non-root branch
 * it inserts a zero-length WGD node, optimizes only q, and measures the AIC
 * improvement.  The branch with the largest improvement above `threshold` is
 * accepted; its WGD node is kept permanently in the tree.  The process repeats
 * until no branch clears the threshold.
 *
 * Re-optimization of base model parameters after accepting a WGD is the
 * caller's responsibility.
 */
class WGDManager {
public:
    struct WGDResult {
        uint   childNodeId;
        double q;
        double lrt;
        double deltaAIC;
    };

public:
    /**
     * @param m        Model parameters.
     * @param tree     The tree (will be modified in-place as WGDs are accepted).
     * @param baseLik  Already-optimized likelihood; parameter values are read
     *                 from here to seed each candidate evaluation.
     * @param threshold Minimum LRT (or ΔAIC) score required to accept a WGD placement.
     *                  Default is 3.841, the χ²(1) critical value at α=0.05.
     */
    WGDManager(ModelParameters* m, std::shared_ptr<PhyloTree> tree, SingleProcessPhyloLikelihood* baseLik, double threshold = 3.841);

    ~WGDManager() {
        for (auto* p : ownedLiks_) {
            LikelihoodUtils::deleteLikelihoodProcess(p);
        }
    }

    void forwardPass();
    const std::vector<WGDResult>& getResults() const { return results_; }
    void printResults() const;
    void writeTree(const std::string& outputPath = "wgd_tree.nwk") const;

private:
    /** Optimize only the WGD.q parameter via Brent. */
    void optimizeQ(SingleProcessPhyloLikelihood* lik);

    /** Read current rate parameter values from a likelihood object. */
    std::map<int, std::vector<double>> extractRateParams(SingleProcessPhyloLikelihood* lik) const;

    /** Clone m_->rDist_ and sync its parameters from the optimized likelihood. */
    std::shared_ptr<DiscreteDistributionInterface> extractRDist(SingleProcessPhyloLikelihood* lik) const;

    // Return relevant node IDs for WGD testing
    std::vector<uint> getCandidates() const;

private:
    ModelParameters* m_;
    std::shared_ptr<PhyloTree> tree_;
    SingleProcessPhyloLikelihood* baseLik_;   // not owned
    double threshold_;

    std::vector<WGDResult> results_;

    // Monotonically increasing edge index — never reuses values so ghost
    // entries left behind by removeSon() never cause setEdgeIndex conflicts.
    uint nextEdgeIdx_ = 0;
    uint nextNodeIdx_ = 0;

    // Likelihood objects created internally after accepting WGDs (owned by us)
    std::vector<SingleProcessPhyloLikelihood*> ownedLiks_;
};

} // namespace bpp

#endif // GENEVOL_WGD_MANAGER_H
