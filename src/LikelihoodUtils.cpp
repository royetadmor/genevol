#include "LikelihoodUtils.h"

using namespace bpp;
using namespace std;


SingleProcessPhyloLikelihood* LikelihoodUtils::createLikelihoodProcess(ModelParameters* m, PhyloTree* tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType) {

    // Create substitution process
    int baseNum = -999;
    bool demiOnlyForEven = false;
    auto mapOfNodeIds = LikelihoodUtils::getMapOfNodeIds(tree);
    std::shared_ptr<DiscreteDistribution> gammaDist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree);
    
    // Create substitution model
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(m->alphabet_, rateParams, baseNum, m->countRange_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, demiOnlyForEven);
    subProcesses->addModel(chrModel, mapOfNodeIds[1]);
    SubstitutionProcess* nsubPro = subProcesses->clone();
    
    // Create likelihood object
    Context* context = new Context();
    bool weightedRootFreqs = true;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *m->container_->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;
}


std::map<uint, vector<uint>> LikelihoodUtils::getMapOfNodeIds(PhyloTree* tree) {
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
    auto sequenceData = lik->getData();
    auto process = &(lik->getSubstitutionProcess());
    auto context = &(lik->getContext());
    delete process;
    delete sequenceData;
    delete context;
    delete lik;
}

void LikelihoodUtils::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, std::vector<std::string> namesAllParams, std::string suffix){
    for (size_t i = 0; i < namesAllParams.size(); i++){
        string cleanParamName;
        if (suffix != ""){
            cleanParamName = namesAllParams[i].substr(0, namesAllParams[i].length() - suffix.length());
        }else{
            cleanParamName = namesAllParams[i];
        }
        uint modelId = getModelFromParamName(cleanParamName);
        int type = getTypeOfParamFromParamName(cleanParamName);
        //should get the type
        
        typeWithParamNames[type][modelId].push_back(namesAllParams[i]);
        if (paramNameAndType){
            (*paramNameAndType)[namesAllParams[i]] = std::pair<int, uint>(type, modelId);
        }
    }

}

uint LikelihoodUtils::getModelFromParamName(string name){
    std::regex modelPattern ("_([\\d]+)");
    std::smatch sm;
    std::regex_search(name, sm, modelPattern);
    std::string modelSuffix = sm[sm.size()-1];
    uint modelId = static_cast<uint>(stoi(modelSuffix));
    return modelId;

}

int LikelihoodUtils::getTypeOfParamFromParamName(string name){
    int type;
    std::map<std::string, int> typeGeneralName;
    updateWithTypeAndCorrespondingName(typeGeneralName);
    auto itParamType = typeGeneralName.begin();
    while(itParamType != typeGeneralName.end()){
        string pattern = itParamType->first;
        if (name.find(pattern) != string::npos){
            type = typeGeneralName[pattern];
            break;

        }
        itParamType ++;
    }
    return type;
}

void LikelihoodUtils::updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
}

int LikelihoodUtils::getParamIndex(string name) {
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
    // find underscore
    size_t pos = paramName.find('_');
    if (pos == std::string::npos)
        throw std::runtime_error("No underscore found in paramName: " + paramName);

    // extract substring before underscore
    std::string prefix = paramName.substr(0, pos-1);

    // collect matches
    std::vector<std::string> result;
    for (const auto& s : listOfParams)
    {
        if (s.find(prefix) != std::string::npos)
        {
            result.push_back(s);
        }
    }

    return result;
}

SingleProcessPhyloLikelihood* LikelihoodUtils::createMyLikelihoodProcess(ModelParameters* m, PhyloTree* tree, std::map<int, std::vector<double>> rateParams, std::vector<int> rateChangeType) {

    // Create substitution process
    auto mapOfNodeIds = LikelihoodUtils::getMapOfNodeIds(tree);
    std::shared_ptr<DiscreteDistribution> gammaDist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree);
    
    // Create substitution model
    std::shared_ptr<GeneCountSubstitutionModel> chrModel = std::make_shared<GeneCountSubstitutionModel>(m->alphabet_, rateParams, GeneCountSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, m);
    subProcesses->addModel(chrModel, mapOfNodeIds[1]);
    SubstitutionProcess* nsubPro = subProcesses->clone();
    
    // Create likelihood object
    Context* context = new Context();
    bool weightedRootFreqs = true;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *m->container_->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;
}

