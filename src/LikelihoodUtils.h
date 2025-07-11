#ifndef GENEVOL_LIKELIHOODUTILS_H
#define GENEVOL_LIKELIHOODUTILS_H

#include <iostream>
#include <string> 
#include <set>

#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>

#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>


#include "ChromosomeSubstitutionModel.h"
#include "ModelParameters.h"



using namespace std;
namespace bpp{
  class LikelihoodUtils{
    public:
      LikelihoodUtils() {}
      virtual ~LikelihoodUtils() {}

    public:
        static SingleProcessPhyloLikelihood* createLikelihoodProcess(ModelParameters* m, PhyloTree* tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType);
        static std::map<uint, vector<uint>> getMapOfNodeIds(PhyloTree* tree);
        static void deleteLikelihoodProcess(SingleProcessPhyloLikelihood* lik);
  };
}

#endif // GENEVOL_LIKELIHOODUTILS_H