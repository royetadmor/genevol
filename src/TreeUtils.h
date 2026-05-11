#ifndef GENEVOL_TREEUTILS_H
#define GENEVOL_TREEUTILS_H

#include <string> 
#include <tuple>
#include <sstream>



#include "ModelParameters.h"

using namespace std;
using namespace bpp;


struct WGDInsertion {
    std::shared_ptr<bpp::PhyloNode> wgdUpper;  // connected to original parent by L/2 edge
    std::shared_ptr<bpp::PhyloNode> wgdLower;  // connected to original child by L/2 edge; zero-length edge above it
    double origLen;
};

class TreeUtils {
    public:
        static double getTreeScalingFactor(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree);
        static int countUniqueStates(const Site site);
        static double getTotalLength(std::shared_ptr<bpp::PhyloTree> tree);
        static void printTopology(std::shared_ptr<bpp::PhyloTree> tree);

        static WGDInsertion insertWGDNode(std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<bpp::PhyloNode> child,
                      uint& nextNodeIdx, uint& nextEdgeIdx);

        static void removeWGDNode(std::shared_ptr<bpp::PhyloTree> tree,
                                  const WGDInsertion& ins, uint& nextEdgeIdx);
};


#endif // GENEVOL_TREEUTILS_H