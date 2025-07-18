#include "MixtureModelLikelihoodFunction.h"


double AlphaLikelihoodFunction::calculateFunctionValue() const {
    double totalMMLikelihood = 0;
    double probabilityPrior = 1 / double(categories_*categories_);
    auto likelihoods = getLikelihoodProcesses();

    for (const auto& lik : likelihoods) {
        double value = lik->getValue() * probabilityPrior;
        totalMMLikelihood += value;
        LikelihoodUtils::deleteLikelihoodProcess(lik);
    }
    return totalMMLikelihood;
}

std::map<std::string, std::vector<double>> AlphaLikelihoodFunction::generateMMValues(int categories, double alphaGain, double betaGain, double alphaLoss, double betaLoss) const {
    std::vector<double> gainValues;
    std::vector<double> lossValues;

    GammaDiscreteDistribution gainGammaDist(categories, alphaGain, betaGain);
    GammaDiscreteDistribution lossGammaDist(categories, alphaLoss, betaLoss);

    for (int i = 0; i < categories; i++)
        gainValues.push_back(gainGammaDist.getCategory(i));

    for (int i = 0; i < categories; i++)
        lossValues.push_back(lossGammaDist.getCategory(i));

    return { {"gain", gainValues}, {"loss", lossValues} };
}

std::vector<SingleProcessPhyloLikelihood*> AlphaLikelihoodFunction::getLikelihoodProcesses() const {
    double alphaGain = getParameterValue("alphaGain0_1");
    double alphaLoss = getParameterValue("alphaLoss0_1");
    double dupl = getParameterValue("dupl0_1");
    double betaGain = getParameterValue("betaGain0_1");
    double betaLoss = getParameterValue("betaLoss0_1");
    std::vector<SingleProcessPhyloLikelihood*> likelihoodProcesses;
    
    // Generate gain/loss values per category
    auto mmValues = generateMMValues(categories_, alphaGain, betaGain, alphaLoss, betaLoss);
    for (double gainValue : mmValues["gain"]) {
        for (double lossValue : mmValues["loss"]) {
            std::map<int, std::vector<double>> MMparamMap = {
                {1, {-999}},               // BaseNumR
                {2, {dupl}},               // Dupl
                {3, {lossValue}},          // Loss
                {4, {gainValue}},          // Gain
                {5, {-999}}                // Demi
            };

            // Build likelihood process
            auto lik = LikelihoodUtils::createLikelihoodProcess(m_, tree_, MMparamMap);
            likelihoodProcesses.push_back(lik);
        }
    }
    return likelihoodProcesses;
}

