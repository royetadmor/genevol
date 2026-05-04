#ifndef GENEVOL_TREEUTILS_H
#define GENEVOL_TREEUTILS_H

#include <string> 
#include <tuple>
#include <sstream>



#include "ModelParameters.h"

using namespace std;
using namespace bpp;


class TreeUtils {
    public:
        static double getTreeScalingFactor(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree);
        static int countUniqueStates(const Site site);
        static double getTotalLength(std::shared_ptr<bpp::PhyloTree> tree);
        static void printTopology(std::shared_ptr<bpp::PhyloTree> tree);

        static std::pair<std::shared_ptr<bpp::PhyloNode>, double>
        insertWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> child,
                      uint& nextNodeIdx, uint& nextEdgeIdx);

        static void removeWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> wgdNode,
                                  double origLen, uint& nextEdgeIdx);
};


#endif // GENEVOL_TREEUTILS_H