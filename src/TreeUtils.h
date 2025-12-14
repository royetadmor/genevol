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
};


#endif // GENEVOL_TREEUTILS_H