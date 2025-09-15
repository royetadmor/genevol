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
#include "ModelParameters.h"
#include "LikelihoodUtils.h"
#include "MixtureModelLikelihoodFunction.h"
#include "GeneCountManager.h"
#include "GeneCountDependencyFunction.h"
#include "GeneCountSubstitutionModel.h"

using namespace bpp;
using namespace std;

void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, ModelParameters* m,double tol, unsigned int maxNumOfIterations, bool mixed=false, unsigned curentIterNum=0);
double getTreeScalingFactor(ModelParameters* m, PhyloTree* tree);
int countUniqueStates(const Site site);

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
    auto rateChangeType = m->rateChangeType_;
    

    // Calculate new likelihood
    auto myNewLik = LikelihoodUtils::createMyLikelihoodProcess(m, tree_, paramMap, rateChangeType);
    std::cout << "New Likelihood: " << myNewLik->getValue() << std::endl;
    if(std::isinf(myNewLik->getValue())) {
        std::cout << "Likelihood is inf, exiting" << std::endl;
        return 1;
    }
    std::cout << "Starting optimization for new model" << std::endl;
    optimizeModelParametersOneDimension(myNewLik, m, 0.1, 2);

    // Create mixture model, calculate likelihood and optimize
    auto geneCountManager = std::make_shared<GeneCountManager>(m, tree_);
    std::cout << "MM Likelihood: " << geneCountManager->getLikelihood() << std::endl;
    geneCountManager->optimizeMixtureModelParametersOneDimension(0.1, 2);
    std::cout << geneCountManager->getLikelihood() << std::endl;
    std::cout << geneCountManager->calculateAIC() << std::endl;

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
    std::vector<int> rateChangeType = m->rateChangeType_;


    vector<string> parametersNames = likelihoodProcess->getSubstitutionModelParameters().getParameterNames();
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
            std::cout << "Previous value of "+ nameOfParam + " is: "+ std::to_string(params.getParameter(nameOfParam).getValue()) << std::endl;

            // This checks if there's a param we don't need to optimize (==fixed param)
            // if (std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), rateParamType)){
            //     continue;
            // }

            // param names corresponding to the parameter type
            std::vector<string> paramsNames = LikelihoodUtils::filterParamsByName(parametersNames, nameOfParam);

            size_t index = LikelihoodUtils::getParamIndex(nameOfParam);
            GeneCountDependencyFunction::FunctionType funcType = static_cast<GeneCountDependencyFunction::FunctionType>(rateChangeType[GeneCountSubstitutionModel::getParamIndexByName(nameOfParam)]);
            GeneCountDependencyFunction* functionOp;
            functionOp = compositeParameter::getDependencyFunction(funcType);

            functionOp->setDomainsIfNeeded(minDomain, maxDomain);
            functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, maxDomain);
            functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);


            delete functionOp;
            std::cout << "Parameter name is: " + nameOfParam << std::endl;
            optimizer->getStopCondition()->setTolerance(tol);
            optimizer->setInitialInterval(lowerBound, upperBound);            
            optimizer->init(params.createSubList(nameOfParam));
            currentLikelihood = optimizer->optimize();
            std::cout << "\nCurrent likelihood: " + std::to_string(currentLikelihood) << std::endl;
            std::cout << nameOfParam + " parameter value after optimization "+ std::to_string(likelihoodProcess->getParameters().getParameter(nameOfParam).getValue()) << std::endl;
        }

        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }

    }
    delete optimizer;
}
