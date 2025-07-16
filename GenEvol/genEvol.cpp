#include <iostream>
#include <string> 
#include <set>

// From bpp-core
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/App/BppApplication.h>

// From bpp-phyl
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>

// From bpp-seq
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// Local modules
#include "ChromosomeSubstitutionModel.h"
#include "ModelParameters.h"
#include "LikelihoodUtils.h"
#include "MixtureModelLikelihoodFunction.h"
#include "GeneCountMixtureModel.h"

using namespace bpp;
using namespace std;

void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, ModelParameters* m,double tol, unsigned int maxNumOfIterations, bool mixed=false, unsigned curentIterNum=0);
double getTreeScalingFactor(ModelParameters* m, PhyloTree* tree);
int countUniqueStates(const Site site);
// void optimizeMixtureModelParametersOneDimension(AlphaLikelihoodFunction* f, ModelParameters* m,double tol, unsigned int maxNumOfIterations, bool mixed=false, unsigned curentIterNum=0);

int main(int args, char **argv) {
    // Set model data and parameters
    BppApplication GenEvol(args, argv, "GenEvol");
    ModelParameters* m = new ModelParameters(GenEvol);

    // Get tree and rescale it
    Newick reader;
    PhyloTree* tree_ = reader.readPhyloTree(m->treeFilePath_);
    double scale_tree_factor = getTreeScalingFactor(m, tree_);
    tree_->scaleTree(scale_tree_factor);

    // Define substitution parameters
    auto paramMap = m->paramMap_;

    // Calculate normal likelihood
    auto newLik = LikelihoodUtils::createLikelihoodProcess(m, tree_, paramMap);
    std::cout << "Likelihood: " << newLik->getValue() << std::endl;
    std::cout << "Starting optimization" << std::endl;
    optimizeModelParametersOneDimension(newLik, m, 0.1, 2);

    // Create mixture model, calculate likelihood and optimize
    auto gcMM = std::make_shared<GeneCountMixtureModel>(m, tree_);
    std::cout << gcMM->getLikelihood() << std::endl;
    gcMM->optimizeMixtureModelParametersOneDimension(0.1, 2);
    std::cout << gcMM->getLikelihood() << std::endl;
    std::cout << gcMM->getParameterValueByName("alphaGain0_1") << std::endl;

    GenEvol.done();
    return 0;
}

int countUniqueStates(const Site site) {
    std::set<int> uniqueStates;
    for (size_t j = 0; j < site.size(); j++) {
        uniqueStates.insert(site[j]); 
    }
    return uniqueStates.size();
}

void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, ModelParameters* m,double tol, unsigned int maxNumOfIterations, bool mixed, unsigned curentIterNum)
{
    // Initialize optimizer
    DerivableSecondOrder* f = likelihoodProcess;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);

    // Can use BRACKET_INWARD or BRACKET_OUTWARD instead
    optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);

    // initializing the likelihood values
    double currentLikelihood = likelihoodProcess->getValue();
    double prevLikelihood;
    int minDomain = m->minState_;
    int maxDomain = m->maxState_;
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
    std::vector<int> rateChangeType = m->rateChangeType_;


    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model

    vector<string> parametersNames = likelihoodProcess->getSubstitutionModelParameters().getParameterNames();
    for (const auto& s : parametersNames) {
        std::cout << s << std::endl;
    }
    LikelihoodUtils::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, parametersNames, 0, "");
    ParameterList params = likelihoodProcess->getParameters();
    size_t nbParams = parametersNames.size();

    // starting iterations of optimization
    for (size_t i = 0; i < maxNumOfIterations; i++)
    {
        prevLikelihood = currentLikelihood;
        for (size_t j = 0; j < nbParams; j ++)
        {
            double lowerBound, upperBound;            
            const string nameOfParam = parametersNames[j];
            int rateParamType = paramNameAndType[nameOfParam].first;
            std::cout << "Previous value of "+ nameOfParam + " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;

            // This checks if there's a param we don't need to optimize (==fixed param)
            // if (std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), rateParamType)){
            //     continue;
            // }

            // param names corresponding to the parameter type
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            Parameter param = params.getParameter(nameOfParam);
            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(rateChangeType[rateParamType - startCompositeParams]);
            ChromosomeNumberDependencyFunction* functionOp;
            functionOp = compositeParameter::setDependencyFunction(funcType);

            functionOp->setDomainsIfNeeded(minDomain, maxDomain);
            functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, maxDomain);
            functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);

            delete functionOp;
            std::cout << "Parameter name is: " + nameOfParam << std::endl;
            optimizer->getStopCondition()->setTolerance(tol);
            optimizer->setInitialInterval(lowerBound, upperBound);            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout << "\nCurrent likelihood: " + std::to_string(currentLikelihood) << std::endl;
            std::cout << nameOfParam + " parameter value after optimization "+ std::to_string(likelihoodProcess->getParameters().getParameter(param.getName()).getValue()) << std::endl;
        }

        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        optimizer->getNumberOfEvaluations();

    }
    delete optimizer;
}

double getTreeScalingFactor(ModelParameters* m, PhyloTree* tree) {
    VectorSiteContainer* container = m->container_;
    if (m->branchMul_ != -999) {
        return m->branchMul_;
    }
    int uniqueStateCount = 0;
    // Iterate over all sites (columns)
    for (size_t i = 0; i < container->getNumberOfSites(); i++) {
        const Site site = container->getSite(i);
        uniqueStateCount += countUniqueStates(site);
    }
    double avgUniqueStateCount = uniqueStateCount/container->getNumberOfSites();
    return avgUniqueStateCount/tree->getTotalLength();
}