#include "TreeUtils.h"

#include <Bpp/Phyl/Tree/PhyloNode.h>
#include <Bpp/Phyl/Tree/PhyloBranch.h>

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

void TreeUtils::printTopology(std::shared_ptr<bpp::PhyloTree> tree)
{
    std::cout << "Tree topology:" << std::endl;
    for (auto& node : tree->getAllNodes()) {
        if (tree->getNodeIndex(node) == tree->getRootIndex()) continue;
        auto edge = tree->getEdgeToFather(node);
        std::cout << "  edge " << tree->getEdgeIndex(edge)
                  << ": node " << tree->getNodeIndex(tree->getFatherOfNode(node))
                  << " -> node " << tree->getNodeIndex(node);
        if (tree->isLeaf(node) && node->hasName())
            std::cout << " (" << node->getName() << ")";
        std::cout << "  len=" << edge->getLength() << std::endl;
    }
}

WGDInsertion
TreeUtils::insertWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> child,
                         uint& nextNodeIdx, uint& nextEdgeIdx, double t)
{
    if (!child)
        throw std::runtime_error("insertWGDNode: child is null.");

    auto parent = tree->getFatherOfNode(child);
    if (!parent)
        throw std::runtime_error("insertWGDNode: cannot insert above root.");

    auto edgeToParent = tree->getEdgeToFather(child);
    if (!edgeToParent)
        throw std::runtime_error("insertWGDNode: child has no incoming edge.");

    double origLen = edgeToParent->getLength();
    double upper = origLen * t;
    double lower = origLen * (1 - t);

    // Detach child from parent
    tree->removeSon(parent, child);

    // Insert wgdUpper between parent and the WGD site (upper half of original branch)
    auto wgdUpper      = std::make_shared<bpp::PhyloNode>("WGD");
    auto upperBranch   = std::make_shared<bpp::PhyloBranch>(upper);
    tree->createNode(parent, wgdUpper, upperBranch);

    // Insert wgdLower below wgdUpper with a zero-length edge (this is the WGD event edge)
    auto wgdLower      = std::make_shared<bpp::PhyloNode>("WGD");
    auto zeroBranch    = std::make_shared<bpp::PhyloBranch>(0.0);
    tree->createNode(wgdUpper, wgdLower, zeroBranch);

    // Re-attach child below wgdLower (lower half of original branch)
    auto lowerBranch   = std::make_shared<bpp::PhyloBranch>(lower);
    tree->link(wgdLower, child, lowerBranch);

    // Assign indices to any new nodes and edges
    for (auto& n : tree->getAllNodes()) {
        if (!tree->hasNodeIndex(n))
            tree->setNodeIndex(n, nextNodeIdx++);
    }
    for (auto& n : tree->getAllNodes()) {
        if (tree->getNodeIndex(n) == tree->getRootIndex()) continue;
        auto e = tree->getEdgeToFather(n);
        if (!tree->hasEdgeIndex(e))
            tree->setEdgeIndex(e, nextEdgeIdx++);
    }

    uint wgdEdgeIdx = tree->getEdgeIndex(tree->getEdgeToFather(wgdLower));
    return {wgdUpper, wgdLower, origLen, wgdEdgeIdx};
}

void TreeUtils::removeWGDNode(std::shared_ptr<bpp::PhyloTree> tree,
                               const WGDInsertion& ins, uint& nextEdgeIdx)
{
    auto parent = tree->getFatherOfNode(ins.wgdUpper);
    if (!parent)
        throw std::runtime_error("removeWGDNode: wgdUpper has no parent.");

    auto sons = tree->getSons(ins.wgdLower);
    if (sons.size() != 1)
        throw std::runtime_error("removeWGDNode: wgdLower expected exactly one child.");
    auto child = sons[0];

    // Detach child from wgdLower, wgdLower from wgdUpper, wgdUpper from parent
    tree->removeSon(ins.wgdLower, child);
    tree->removeSon(ins.wgdUpper, ins.wgdLower);
    tree->removeSon(parent, ins.wgdUpper);

    // Restore original edge
    auto restoredBranch = std::make_shared<bpp::PhyloBranch>(ins.origLen);
    tree->link(parent, child, restoredBranch);

    tree->deleteNode(ins.wgdUpper);
    tree->deleteNode(ins.wgdLower);

    if (!tree->hasEdgeIndex(restoredBranch))
        tree->setEdgeIndex(restoredBranch, nextEdgeIdx++);
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