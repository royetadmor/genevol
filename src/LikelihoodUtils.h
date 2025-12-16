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
#include <Bpp/Seq/Container/AlignmentData.h>

#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>

#include "ModelParameters.h"
#include "GeneCountSubstitutionModel.h"



using namespace std;
namespace bpp{
  class LikelihoodUtils{
    public:
      LikelihoodUtils() {}
      virtual ~LikelihoodUtils() {}

    public:
        static std::map<uint, vector<uint>> getMapOfNodeIds(std::shared_ptr<bpp::PhyloTree> tree);
        static void deleteLikelihoodProcess(SingleProcessPhyloLikelihood* lik);
        static int getParamIndex(string name);
        static std::vector<string> filterParamsByName(std::vector<std::string> listOfParams, std::string paramName);
        static SingleProcessPhyloLikelihood* createLikelihoodProcess(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType, std::map<string, string> constraintedParams, std::shared_ptr<DiscreteDistributionInterface> rDist);
        static void setProcessConstraintedParams(std::map<string, string> constraintedParams, AbstractParameterAliasable* process);
        static bool isFixedParam(const std::string& name, const std::vector<string> params);
        static std::vector<double> calculateExpectedRatePerSite(SingleProcessPhyloLikelihood* lik, const bool normalize);
    private:
        static vector<string> getParametersByName(ParameterList params, string name);
        static void normalizeVector(vector<double>& data);
  };
}

#endif // GENEVOL_LIKELIHOODUTILS_H