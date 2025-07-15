#include "MixtureModelLikelihoodFunction.h"


double AlphaLikelihoodFunction::calculateFunctionValue() const {
    double alpha_gain = getParameterValue("alphaGain0_1");
    double alpha_loss = getParameterValue("alphaLoss0_1");
    double dupl = getParameterValue("dupl0_1");
    double beta_gain = 1.0 / alpha_gain;
    double probabilityPrior = 1/double(categories_*categories_);
    double totalMMLikelihood = 0;

    auto mmValues = generateMMValues(categories_, alpha_gain, beta_gain, alpha_loss);
    for (double gainValue: mmValues["gain"])
    {
        for (double lossValue: mmValues["loss"]) {
            
            std::map<int, std::vector<double>> MMparamMap = {
                {1, {-999}}, // BaseNumR
                {2, {dupl}}, // Dupl
                {3, {lossValue, 0.1}}, // Loss
                {4, {gainValue, 0.1}}, // Gain
                {5, {-999}} // Demi
            };
            std::vector<int> rateChangeType = m_->rateChangeType_;
            auto MMLikelihood = LikelihoodUtils::createLikelihoodProcess(m_, tree_, MMparamMap);
            totalMMLikelihood += MMLikelihood->getValue() * probabilityPrior;
            LikelihoodUtils::deleteLikelihoodProcess(MMLikelihood);
        }
    }
    return totalMMLikelihood;
}

std::map<std::string, std::vector<double>> AlphaLikelihoodFunction::generateMMValues(int categories, double alphaGain, double betaGain, double shapeLoss) const {
    std::vector<double> gain_values;
    std::vector<double> loss_values;

    GammaDiscreteDistribution gainGammaDist(categories, alphaGain, betaGain);
    GammaDiscreteDistribution lossGammaDist(categories, shapeLoss, shapeLoss);

    for (int i = 0; i < categories; i++)
        gain_values.push_back(gainGammaDist.getCategory(i));

    for (int i = 0; i < categories; i++)
        loss_values.push_back(lossGammaDist.getCategory(i));

    return { {"gain", gain_values}, {"loss", loss_values} };
}
