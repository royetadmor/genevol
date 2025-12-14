#include "TreeUtils.h"

double TreeUtils::getTotalLength(std::shared_ptr<bpp::PhyloTree> tree) {
    double treeLength = 0;
    Vdouble branchLengths = tree->getBranchLengths();
    for (size_t i = 0; i < branchLengths.size(); i++){
        treeLength += branchLengths[i];
    }
    return treeLength;
}


int TreeUtils::countUniqueStates(const Site site) {
    std::set<int> uniqueStates;
    for (size_t j = 0; j < site.size(); j++) {
        uniqueStates.insert(site[j]); 
    }
    return uniqueStates.size();
}

double TreeUtils::getTreeScalingFactor(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree) {
    auto container = m->container_;
    if (m->branchMul_ != -999) {
        return m->branchMul_;
    }
    int uniqueStateCount = 0;
    // Iterate over all sites (columns)
    for (size_t i = 0; i < container->getNumberOfSites(); i++) {
        const Site site = container->site(i);
        uniqueStateCount += countUniqueStates(site);
    }
    double avgUniqueStateCount = uniqueStateCount/container->getNumberOfSites();
    double totalLength = getTotalLength(tree);
    return avgUniqueStateCount/totalLength;
}