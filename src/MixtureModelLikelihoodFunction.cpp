#include "MixtureModelLikelihoodFunction.h"


double MixtureModelLikelihoodFunction::calculateFunctionValue() const {
    double totalMMLikelihood = 0;
    double probabilityPrior = 1 / double(categories_*categories_);
    auto likelihoods = getLikelihoodProcesses();
    size_t numberOfSites = m_->container_->getNumberOfSites();
    for (size_t i = 0; i < numberOfSites; i++) {
        double familyLikelihood = 0;
        for (const auto& lik : likelihoods) {
            familyLikelihood += bpp::ExtendedFloat::convert(lik->getLikelihoodPerSite()[i]) * probabilityPrior;
        }
        totalMMLikelihood += log(familyLikelihood);
    }

    for (const auto& lik : likelihoods) {
        LikelihoodUtils::deleteLikelihoodProcess(lik);
    }
    return -totalMMLikelihood;
}

std::map<std::string, std::vector<double>> MixtureModelLikelihoodFunction::generateMMValues(int categories, double alphaGain, double betaGain, double alphaLoss, double betaLoss) const {
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

std::vector<SingleProcessPhyloLikelihood*> MixtureModelLikelihoodFunction::getLikelihoodProcesses() const {
    std::vector<SingleProcessPhyloLikelihood*> likelihoodProcesses;
    const auto alphaGain = getParameterValue("alphaGain0_1");
    const auto betaGain = getParameterValue("betaGain0_1");
    const auto alphaLoss = getParameterValue("alphaLoss0_1");
    const auto betaLoss = getParameterValue("betaLoss0_1");
    const auto innovation = getParameterValue("innovation0_1");
    const auto elimination = getParameterValue("elimination0_1");
    const auto& lossRateParams = getRateParametersByEvent(GeneCountSubstitutionModel::LOSS);
    const auto& gainRateParams = getRateParametersByEvent(GeneCountSubstitutionModel::GAIN);
    const auto& innovationRateParams = getRateParametersByEvent(GeneCountSubstitutionModel::INNOVATION);
    const auto& eliminationRateParams = getRateParametersByEvent(GeneCountSubstitutionModel::ELIMINATION);

    // Generate gain/loss values per category
    const auto mmValues = generateMMValues(categories_, alphaGain, betaGain, alphaLoss, betaLoss);

    for (double gainValue : mmValues.at("gain")) {
        for (double lossValue : mmValues.at("loss")) {
            std::vector<double> lossValues{lossValue};
            std::vector<double> gainValues{gainValue};
            std::vector<double> innovationValues{innovation};
            std::vector<double> eliminationValues{elimination};
            lossValues.insert(lossValues.end(), lossRateParams.begin(), lossRateParams.end());
            gainValues.insert(gainValues.end(), gainRateParams.begin(), gainRateParams.end());
            innovationValues.insert(innovationValues.end(), innovationRateParams.begin(), innovationRateParams.end());
            eliminationValues.insert(eliminationValues.end(), eliminationRateParams.begin(), eliminationRateParams.end());

            std::map<int, std::vector<double>> MMparamMap = {
                {GeneCountSubstitutionModel::LOSS, lossValues},
                {GeneCountSubstitutionModel::GAIN, gainValues},
                {GeneCountSubstitutionModel::INNOVATION, innovationValues},
                {GeneCountSubstitutionModel::ELIMINATION, eliminationValues},
            };

            likelihoodProcesses.push_back(
                LikelihoodUtils::createLikelihoodProcess(m_, tree_, MMparamMap, m_->mixtureRateChangeType_, {})
            );
        }
    }

    return likelihoodProcesses;
}


std::vector<double> MixtureModelLikelihoodFunction::getRateParametersByEvent(int eventType) const {
    std::vector<double> values;
    if (rateParams.find(eventType) != rateParams.end()) {
        std::vector<string> paramNames = rateParams.at(eventType);
        for (const auto& name : paramNames) {
            values.push_back(getParameterValue(name));
        }
    }
    return values;
}

