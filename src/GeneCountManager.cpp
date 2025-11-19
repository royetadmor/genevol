#include "GeneCountManager.h"
#include <iostream>

using namespace bpp;
using namespace std;

void GeneCountManager::optimizeMixtureModelParametersOneDimension(double tol, unsigned int maxNumOfIterations, bool mixed, unsigned curentIterNum)
{
    // Initialize optimizer
    auto f = likelihoodFunction_.get();
    BrentOneDimension* optimizer = new BrentOneDimension(likelihoodFunction_);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    optimizer->getStopCondition()->setTolerance(tol);
    // Can use BRACKET_INWARD or BRACKET_OUTWARD instead
    optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    // initializing the likelihood values
    double currentLikelihood = f->getValue();
    double prevLikelihood;
    int minDomain = m_->alphabet_->getMin();
    int maxDomain = m_->alphabet_->getMax();
    std::vector<int> rateChangeType = m_->mixtureRateChangeType_;


    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model

    vector<string> parametersNames = f->getIndependentParameters().getParameterNames();
    ParameterList params = f->getIndependentParameters();
    size_t nbParams = parametersNames.size();

    
    // starting iterations of optimization
    for (size_t i = 0; i < maxNumOfIterations; i++)
    {
        prevLikelihood = currentLikelihood;
        for (size_t j = 0; j < nbParams; j ++)
        {
            double lowerBound, upperBound;            
            const string nameOfParam = parametersNames[j];
            std::cout << "Previous value of "+ nameOfParam + " is: "+ std::to_string(params.getParameter(nameOfParam)->getValue()) << std::endl;

            if (LikelihoodUtils::isFixedParam(nameOfParam, m_->mixtureFixedParams_)) {
                std::cout << "Skipping " << nameOfParam << std::endl;
                continue;
            }
            std::vector<string> paramsNames = LikelihoodUtils::filterParamsByName(parametersNames, nameOfParam);
            size_t index = LikelihoodUtils::getParamIndex(nameOfParam);

            // TODO: this is a hack you should fix. The problem it currently fix is:
            // - I assume every parameter is constant (which is false because I support dependency funcs in MM)
            // - So I check when the index is 1 and set it manually to linear - PROBLEM! need to get the dep. function somehow
            // - since the name of the rate parameter is "rategain1_1" it doesn't associate it with the other gain param (alphaGain)
            // - so the list is missing one element, so when `updateBounds` is called it tries to get the 2nd parameter and fails
            // - because he isn't there. The manual fix is to add another random parameter name from the back so it can find
            // - the right parameter name at the right index, but this is bad and should fix!
            // so to summarize: need to fix the dependecny function and the paramsNames list.
            GeneCountDependencyFunction::FunctionType funcType = GeneCountDependencyFunction::FunctionType::CONSTANT;
            if(index == 1) {
                funcType = GeneCountDependencyFunction::FunctionType::LINEAR;
                paramsNames.insert(paramsNames.begin(), "alphaLoss0_1");
            }
            GeneCountDependencyFunction* functionOp;
            functionOp = compositeParameter::getDependencyFunction(funcType);

            functionOp->setDomainsIfNeeded(minDomain, maxDomain);
            functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, maxDomain);
            functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
            
            delete functionOp;
            std::cout << "Parameter name is: " + nameOfParam << std::endl;
            optimizer->setInitialInterval(lowerBound, upperBound);            
            optimizer->init(params.createSubList(nameOfParam));
            currentLikelihood = optimizer->optimize();
            std::cout << "\nCurrent likelihood: " + std::to_string(currentLikelihood) << std::endl;
            std::cout << nameOfParam + " parameter value after optimization "+ std::to_string(f->getParameters().getParameter(nameOfParam)->getValue()) << std::endl;
        }

        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        optimizer->getNumberOfEvaluations();
    }
    delete optimizer;
}

double GeneCountManager::calculateAIC() {
    auto numOfParams = likelihoodFunction_->getIndependentParameters().size();
    double AIC = 2*(likelihoodFunction_->getValue()) + (2*numOfParams);
    return AIC;
}
