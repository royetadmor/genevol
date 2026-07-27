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
#include <Bpp/Phyl/Likelihood/DataFlow/ForwardLikelihoodTree.h>

#include "ModelParameters.h"
#include "GeneCountSubstitutionModel.h"
#include "NegBinomialFrequencySet.h"
#include "PoissonFrequencySet.h"
#include "WGDSubstitutionModel.h"
#include "ExtendedBrentOptimizer.h"
#include "GeneCountDependencyFunction.h"
#include "ModelAdequacyUtils.h"
#include <Bpp/Numeric/Random/RandomTools.h>



using namespace std;
namespace bpp{
  class LikelihoodUtils{
    public:
      LikelihoodUtils() {}
      virtual ~LikelihoodUtils() {}

    public:
        static void deleteLikelihoodProcess(SingleProcessPhyloLikelihood* lik);
        static int getParamIndex(string name);
        static std::vector<string> filterParamsByName(std::vector<std::string> listOfParams, std::string paramName);
        static std::shared_ptr<NonHomogeneousSubstitutionProcess> createSubstitutionProcess(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType, std::shared_ptr<DiscreteDistributionInterface> rDist, std::map<uint, double> wgdQMap = {}, double qInit = 0.0, double rootLambda = -1.0);
        static SingleProcessPhyloLikelihood* createLikelihoodProcess(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType, std::map<string, string> constraintedParams, std::shared_ptr<DiscreteDistributionInterface> rDist, std::map<uint, double> wgdQMap = {}, double qInit = 0.0);
        static void setProcessConstraintedParams(std::map<string, string> constraintedParams, AbstractParameterAliasable* process);
        static bool isFixedParam(const std::string& name, const std::vector<string> params);
        static std::vector<double> calculateExpectedRatePerSite(SingleProcessPhyloLikelihood* lik, const bool normalize);
        static void printResults(SingleProcessPhyloLikelihood* lik, bool printRate4Site = false);
        static void printRootFreqsPerSite(SingleProcessPhyloLikelihood* lik);
        static void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, ModelParameters* m,double tol, unsigned int maxNumOfIterations);
        static SingleProcessPhyloLikelihood* multiStartOptimize(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::vector<int> rateChangeType, std::map<string, string> constraintedParams, std::shared_ptr<DiscreteDistributionInterface> rDist);
    private:
        static vector<string> getParametersByName(ParameterList params, string name);
        static void normalizeVector(vector<double>& data);
        static std::shared_ptr<PoissonFrequencySet> poissonRootFreqSet(ModelParameters* m, std::shared_ptr<const StateMapInterface> stateMap, double rootLambda = -1.0);
        static std::shared_ptr<NegBinomialFrequencySet> negBinRootFreqSet(ModelParameters* m, std::shared_ptr<const StateMapInterface> stateMap);
        static std::map<int, std::vector<double>> createRandomRateParams(ModelParameters* m);
        static std::vector<SingleProcessPhyloLikelihood*> selectTopK(const std::vector<SingleProcessPhyloLikelihood*>& candidates, int k);
  };
}

#endif // GENEVOL_LIKELIHOODUTILS_H