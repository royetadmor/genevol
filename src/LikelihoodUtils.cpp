#include "LikelihoodUtils.h"

using namespace bpp;
using namespace std;


SingleProcessPhyloLikelihood* LikelihoodUtils::createLikelihoodProcess(ModelParameters* m, std::shared_ptr<bpp::PhyloTree> tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType, std::map<string, vector<string>> constraintedParams, DiscreteDistribution* rDist) {
  
    // Create substitution process components
    std::shared_ptr<DiscreteDistributionInterface> gammaDist = std::shared_ptr<DiscreteDistributionInterface>(rDist);
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree);

    // Create substitution model
    std::shared_ptr<GeneCountSubstitutionModel> subModel = std::make_shared<GeneCountSubstitutionModel>(m->alphabet_, rateParams, GeneCountSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, m);

    // Calculate and set root frequencies
    auto alphabetSize = m->alphabet_->getSize();
    std::vector<double> freq(alphabetSize, (1.0/alphabetSize));
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(subModel->getStateMap(), freq);
    std::shared_ptr<FrequencySetInterface> rootFrequencies = static_pointer_cast<FrequencySetInterface>(rootFreqsFixed);

    // Create substitution process
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree, rootFrequencies);
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

void LikelihoodUtils::setProcessConstraintedParams(std::map<string, vector<string>> constraintedParams, AbstractParameterAliasable* aliasable) {
    ParameterList params = aliasable->getParameters();
    for (const auto& pair : constraintedParams) {
        string p1 = LikelihoodUtils::getParameterByName(params, pair.first);
        for (const string& s : pair.second) {
            string p2 = LikelihoodUtils::getParameterByName(params, s);
            aliasable->aliasParameters(p1,p2);
        }
    }
}

// Given a parameter name, returns the matching parameter name from the substitution model
// i.e., if the parameter name is "gain", it will return "GeneCount.gain0_1" etc.
string LikelihoodUtils::getParameterByName(ParameterList params, string name) {
    for (size_t i = 0; i < params.size(); ++i)
    {
        const string paramName = params[i].getName();
        if (paramName.find(name) != std::string::npos) // substring found
        {
            return paramName;
        }
    }
    const string err = "Parameter name not found: " + name;
    throw std::runtime_error(err);
}

bool LikelihoodUtils::isFixedParam(const std::string& name, const std::vector<string> params) {
    for (const std::string& element : params) {
        if (name.find(element) != std::string::npos) {
            return true;
        }
    }
    return false;
}

