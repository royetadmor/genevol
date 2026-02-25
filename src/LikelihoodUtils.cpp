#include "LikelihoodUtils.h"

using namespace bpp;
using namespace std;


SingleProcessPhyloLikelihood* LikelihoodUtils::createLikelihoodProcess(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType, std::map<string, string> constraintedParams, std::shared_ptr<DiscreteDistributionInterface> rDist) {
    // Create substitution process components
    auto parTree = std::make_shared<ParametrizablePhyloTree>(*tree);

    // Create substitution model
    std::shared_ptr<GeneCountSubstitutionModel> subModel = std::make_shared<GeneCountSubstitutionModel>(m->alphabet_, rateParams, GeneCountSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, m);

    // Calculate and set root frequencies
    // Can switch to getRootFrequncies(m, tree, rDist, subModel);
    auto freq = poissonRootFreq(m); 
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(subModel->getStateMap(), freq);
    std::shared_ptr<FrequencySetInterface> rootFrequencies = static_pointer_cast<FrequencySetInterface>(rootFreqsFixed);

    // Create substitution process
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(rDist, parTree, rootFrequencies);
    auto mapOfNodeIds = LikelihoodUtils::getMapOfNodeIds(tree);
    subProcesses->addModel(subModel, mapOfNodeIds[1]);
    auto nsubPro = std::shared_ptr<bpp::SubstitutionProcessInterface>(subProcesses->clone());

    // Check for constrainted params
    if (!constraintedParams.empty()) {
        auto aliasable = dynamic_cast<AbstractParameterAliasable*>(nsubPro.get());
        LikelihoodUtils::setProcessConstraintedParams(constraintedParams, aliasable);
    }

    // Create likelihood object
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, m->container_, nsubPro);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;
}

std::map<uint, vector<uint>> LikelihoodUtils::getMapOfNodeIds(std::shared_ptr<bpp::PhyloTree> tree) {
    std::map<uint, vector<uint>> mapModelNodesIds;
    auto nodes = tree->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        uint nodeId = tree->getNodeIndex(nodes[i]);
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        mapModelNodesIds[1].push_back(nodeId);
    }
    // mapModelNodesIds[1].push_back(tree->getRootIndex());
    return mapModelNodesIds;
}

void LikelihoodUtils::deleteLikelihoodProcess(SingleProcessPhyloLikelihood* lik) {
    auto context = &(lik->context());
    delete context;
    delete lik;
}

int LikelihoodUtils::getParamIndex(string name) {
    // Special case for gamma distribution parameter names which does not contain an underscore
    if (name.find("Gamma") != std::string::npos) {
        return 0;
    }
    size_t pos = name.find('_');
    if (pos == std::string::npos || pos == 0) {
        throw std::runtime_error("No underscore found in string: " + name);
    }
    char index = name[pos - 1];
    if (!std::isdigit(static_cast<unsigned char>(index))) {
        throw std::runtime_error(std::string("Character before underscore is not a digit: ") + index);
    }

    return index - '0'; // convert char digit to int
}

std::vector<string> LikelihoodUtils::filterParamsByName(std::vector<std::string> listOfParams, std::string paramName)
{
    std::vector<std::string> result;

    // Special case for gamma distribution parameter names which does not contain an underscore
    if (paramName.find("Gamma") != std::string::npos) {
        result.push_back(paramName);
        return result;
    }
    
    // find underscore
    size_t pos = paramName.find('_');
    if (pos == std::string::npos) {
        throw std::runtime_error("No underscore found in paramName: " + paramName);
    }

    // extract substring before underscore
    std::string prefix = paramName.substr(0, pos-1);

    // collect matches
    for (const auto& s : listOfParams)
    {
        if (s.find(prefix) != std::string::npos)
        {
            result.push_back(s);
        }
    }

    return result;
}

void LikelihoodUtils::setProcessConstraintedParams(std::map<string, string> constraintedParams, AbstractParameterAliasable* aliasable) {
    ParameterList params = aliasable->getParameters();
    for (const auto& pair : constraintedParams) {
        vector<string> p1 = LikelihoodUtils::getParametersByName(params, pair.first);
        vector<string> p2 = LikelihoodUtils::getParametersByName(params, pair.second);
        for (const auto& s1 : p1) {
            for (const auto& s2 : p2) {
                if (s1.substr(s1.size() - 3) == s2.substr(s2.size() - 3)) {
                    std::cout << "Aliasing " << s1 << " <-> " << s2 << std::endl;
                    aliasable->aliasParameters(s1,s2);
                    continue;
                }
            }
        }
    }
}

// Given a parameter name, returns all matching parameters names from the substitution model
// i.e., if the parameter name is "gain", it will return "GeneCount.gain0_1, GeneCount.gain1_1" etc.
vector<string> LikelihoodUtils::getParametersByName(ParameterList params, string name) {
    vector<string> res;
    for (size_t i = 0; i < params.size(); ++i)
    {
        const string paramName = params[i].getName();
        if (paramName.find(name) != std::string::npos) // substring found
        {
            res.push_back(paramName);
        }
    }
    return res;
}

std::vector<double> LikelihoodUtils::calculateExpectedRatePerSite(SingleProcessPhyloLikelihood* lik, const bool normalize) {
  std::vector<double> res;
  std::shared_ptr<const SubstitutionProcessInterface> pSP = lik->getSubstitutionProcess();
  auto pDD = pSP->getRateDistribution();
  auto sites = lik->getData();
  for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
    auto post = lik->getPosteriorProbabilitiesForSitePerClass(i);
    double expectedR = 0;
    for (size_t j = 0; j < post.size(); ++j) {
      auto cat = pDD->getCategory(j);
      expectedR += post[j] * cat;
    }
    res.push_back(expectedR);
  }
  if (normalize) {
    normalizeVector(res);
  }
  return res;
}

void LikelihoodUtils::normalizeVector(vector<double>& data) {
    if (data.empty()) {
        throw invalid_argument("Input vector is empty.");
    }

    // Compute mean
    double sum = accumulate(data.begin(), data.end(), 0.0);
    double size = static_cast<double>(data.size());
    double mean = sum / size;

    // Compute standard deviation
    double sq_sum = 0.0;
    for (double val : data) {
        sq_sum += (val - mean) * (val - mean);
    }
    double stdev = sqrt(sq_sum / size);

    if (stdev == 0.0) {
        throw runtime_error("Standard deviation is zero. Cannot normalize.");
    }

    // Normalize each element
    for (double& val : data) {
        val = (val - mean) / stdev;
    }
}

bool LikelihoodUtils::isFixedParam(const std::string& name, const std::vector<string> params) {
    for (const std::string& element : params) {
        if (name.find(element) != std::string::npos) {
            return true;
        }
    }
    return false;
}

double LikelihoodUtils::calculateAIC(SingleProcessPhyloLikelihood* lik) {
    auto numOfParams = lik->getSubstitutionProcess()->getIndependentParameters().size();
    double AIC = 2*(lik->getValue()) + (2*numOfParams);
    return AIC;
}

void LikelihoodUtils::printResults(SingleProcessPhyloLikelihood* lik, bool printRate4Site) {
    double likelihood = lik->getValue();
    double aic = LikelihoodUtils::calculateAIC(lik);

    std::cout << "Likelihood: " << likelihood << std::endl;
    std::cout << "AIC score: " << aic << std::endl;

    ParameterList params = lik->getParameters();
    std::cout << "Parameter values:" << std::endl;
    for (size_t i = 0; i < params.size(); ++i) {
        if (params[i].getName().find("BrLen") == std::string::npos) { // Print non-branch length params
            std::cout << "  " << params[i].getName() << " = " << params[i].getValue() << std::endl;
        }
    }

    if (printRate4Site) {
        std::cout << "Calculating expected rate per site" << std::endl;
        auto ratePerSite = LikelihoodUtils::calculateExpectedRatePerSite(lik, false);
        for (size_t i = 0; i < ratePerSite.size(); ++i) {
            std::cout << "Expected rate at site " << i << " is " << ratePerSite[i] << std::endl;
        }
    }
}

std::vector<double> LikelihoodUtils::getRootFrequncies(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::shared_ptr<DiscreteDistributionInterface> rDist, std::shared_ptr<GeneCountSubstitutionModel> model) {
    size_t nbSite = m->container_->getNumberOfSites();
    size_t nbState = m->alphabet_->getSize();
    size_t nbCat = rDist->getNumberOfCategories();
    std::vector<double> res(nbState, 0.0);
    auto parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    Context* context = new Context();

    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(rDist, parTree);
    auto mapOfNodeIds = LikelihoodUtils::getMapOfNodeIds(tree);
    subProcesses->addModel(model, mapOfNodeIds[1]);
    auto nsubPro = std::shared_ptr<bpp::SubstitutionProcessInterface>(subProcesses->clone());

    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, m->container_, nsubPro);
    lik->makeLikelihoods();

    for (size_t cat = 0; cat < nbCat; ++cat) {
        double prior = rDist->getProbability(cat);
        auto flt = lik->getForwardLikelihoodTree(cat);
        auto sumOfWeights = CWiseAdd<RowLik, MatrixLik>::create(*context, {flt->getForwardLikelihoodArrayAtRoot()}, RowVectorDimension(Eigen::Index(nbSite)));
        auto rootSumOfWeights = CWiseAdd<DataLik, RowLik>::create(*context, {sumOfWeights}, Dimension<DataLik>());
        auto converted = Convert<MatrixLik, Transposed<MatrixLik>>::create(*context, {flt->getForwardLikelihoodArrayAtRoot()}, MatrixDimension ((size_t)nbSite, nbState));
        auto rootCondLik = CWiseAdd<RowLik, MatrixLik>::create(*context, {converted}, RowVectorDimension(Eigen::Index(nbState)));
        auto freqs_ef = CWiseDiv<RowLik, std::tuple<RowLik, DataLik>>::create(*context,{rootCondLik, rootSumOfWeights}, RowVectorDimension(Eigen::Index(nbState)));
        auto rFreqs = Convert<Eigen::RowVectorXd, ExtendedFloatRowVectorXd>::create(*context, {freqs_ef}, RowVectorDimension (Eigen::Index (nbState)));
        for (size_t i = 0; i < nbState; ++i)
        {
            res[i] += prior * rFreqs->targetValue()(i);
        }
    }
    return res;
}

// Function mostly for debugging, consider removing in the future
void LikelihoodUtils::printRootFreqsPerSite(SingleProcessPhyloLikelihood* lik)
{
    auto process = lik->getSubstitutionProcess();

    if (process->hasRootFrequencySet())
    {
    auto rootFreqSet = process->getRootFrequencySet();
    const std::vector<double>& freqs = rootFreqSet->getFrequencies();

    std::cout << "Root frequencies:\n";
    for (size_t i = 0; i < freqs.size(); ++i)
        std::cout << "  f[" << i << "] = " << freqs[i] << '\n';
    }
}

std::vector<double> LikelihoodUtils::poissonRootFreq(ModelParameters* m) {
    double total = 0.0;
    double totalDataSize = 0.0;
    size_t nbState = m->alphabet_->getSize();
    std::vector<double> rFreq(nbState, 0.0);
    const size_t nSites = m->container_->getNumberOfSites();

    for (size_t i = 0; i < nSites; ++i)
    {
        const bpp::Site& site = m->container_->site(i);
        const size_t siteSize = site.size();
        totalDataSize += siteSize;
        for (size_t j = 0; j < siteSize; ++j)
        {
            total += site[j];
        }
    }
    double lambda = total/totalDataSize;
    auto pDist = std::make_shared<PoissonDistribution>(lambda, nbState);
    for (size_t i = 0; i < nbState; i++)
    {
        rFreq[i] = pDist->probability(i);
    }

    return rFreq;
}