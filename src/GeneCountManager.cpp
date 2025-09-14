#include "GeneCountManager.h"
#include <iostream>

using namespace bpp;
using namespace std;

void GeneCountManager::optimizeMixtureModelParametersOneDimension(double tol, unsigned int maxNumOfIterations, bool mixed, unsigned curentIterNum)
{
    // Initialize optimizer
    auto f = likelihoodFunction_.get();
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    optimizer->getStopCondition()->setTolerance(tol);
    // Can use BRACKET_INWARD or BRACKET_OUTWARD instead
    optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);

    // initializing the likelihood values
    double currentLikelihood = f->getValue();
    double prevLikelihood;
    int minDomain = m_->alphabet_->getMin();
    int maxDomain = m_->alphabet_->getMax();
    std::vector<int> rateChangeType = m_->mixtureRateChangeType_;


    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model

    vector<string> parametersNames = f->getParameters().getParameterNames();
    ParameterList params = f->getParameters();
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
            std::vector<string> paramsNames = LikelihoodUtils::filterParamsByName(parametersNames, nameOfParam);
            // Parameter param = params.getParameter(nameOfParam);

            size_t index = LikelihoodUtils::getParamIndex(nameOfParam);
            // Since we don't estimate the parameters directly from the data, but rather from the gamma distrbution,
            // the parameters value's don't depend on the current state, hence constant.
            GeneCountDependencyFunction::FunctionType funcType = GeneCountDependencyFunction::FunctionType::CONSTANT;
            GeneCountDependencyFunction* functionOp;
            functionOp = NcompositeParameter::getDependencyFunction(funcType);

            functionOp->setDomainsIfNeeded(minDomain, maxDomain);
            functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, maxDomain);
            functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);

            delete functionOp;
            std::cout << "Parameter name is: " + nameOfParam << std::endl;
            optimizer->setInitialInterval(lowerBound, upperBound);            
            optimizer->init(params.createSubList(nameOfParam));
            currentLikelihood = optimizer->optimize();
            std::cout << "\nCurrent likelihood: " + std::to_string(currentLikelihood) << std::endl;
            std::cout << nameOfParam + " parameter value after optimization "+ std::to_string(f->getParameters().getParameter(nameOfParam).getValue()) << std::endl;
        }

        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        optimizer->getNumberOfEvaluations();

    }
    delete optimizer;
}

double GeneCountManager::calculateAIC() {
    auto numOfParams = likelihoodFunction_->getParametersCount();
    double AIC = 2*(likelihoodFunction_->getValue()) + (2*numOfParams);
    return AIC;
}
