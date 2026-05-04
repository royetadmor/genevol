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

std::pair<std::shared_ptr<bpp::PhyloNode>, double>
TreeUtils::insertWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> child,
                         uint& nextNodeIdx, uint& nextEdgeIdx)
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

    tree->removeSon(parent, child);

    auto wgdNode    = std::make_shared<bpp::PhyloNode>("WGD");
    auto zeroBranch = std::make_shared<bpp::PhyloBranch>(0.0);
    tree->createNode(parent, wgdNode, zeroBranch);

    auto newBranch = std::make_shared<bpp::PhyloBranch>(origLen);
    tree->link(wgdNode, child, newBranch);

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

    return {wgdNode, origLen};
}

void TreeUtils::removeWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> wgdNode,
                               double origLen, uint& nextEdgeIdx)
{
    auto parent = tree->getFatherOfNode(wgdNode);
    if (!parent)
        throw std::runtime_error("removeWGDNode: WGD node has no parent.");

    auto sons = tree->getSons(wgdNode);
    if (sons.size() != 1)
        throw std::runtime_error("removeWGDNode: expected exactly one child.");
    auto child = sons[0];

    tree->removeSon(wgdNode, child);
    tree->removeSon(parent, wgdNode);

    auto restoredBranch = std::make_shared<bpp::PhyloBranch>(origLen);
    tree->link(parent, child, restoredBranch);

    tree->deleteNode(wgdNode);

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