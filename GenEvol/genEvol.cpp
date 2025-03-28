#include <iostream>
#include <string> 
#include <set>

#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/AutoParameter.h>

#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>


#include "ChromosomeSubstitutionModel.h"


using namespace bpp;
using namespace std;


std::tuple<int, int> getAlphabetLimit(const std::string& filePath);
int getGeneCountRange(const VectorSiteContainer* container);
VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
std::vector<std::string> getSpeciesListFromFile(const std::string& filePath);
std::map<uint, vector<uint>> getMapOfNodeIds(PhyloTree* tree);
void updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, std::vector<std::string> namesAllParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::string suffix);
void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, double tol, unsigned int maxNumOfIterations, bool mixed=false, unsigned curentIterNum=0);
uint getModelFromParamName(string name);
void updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName);
int getTypeOfParamFromParamName(string name);
double getTreeScalingFactor(const VectorSiteContainer* container, PhyloTree* tree);
int countUniqueStates(const Site* site);

int main(int args, char **argv) {
    // Get tree path and count path from args
    std::string treePath = "test_data/tree.newick";
    std::string dataPath = "test_data/data.fasta";

    // Set Alphabet
    int minState, maxState;
    tie(minState, maxState) = getAlphabetLimit(dataPath);
    IntegerAlphabet* alphabet = new IntegerAlphabet(110, 1); // TODO: hardcoded

    // Get sequences
    VectorSiteContainer* container = readGeneFamilyFile(dataPath, alphabet);
    std::vector<std::string> species = getSpeciesListFromFile(dataPath);

    // Get tree and rescale it
    Newick reader;
    PhyloTree* tree_ = reader.readPhyloTree(treePath);
    int countRange = getGeneCountRange(container);
    countRange = 50; // TODO: hardcoded
    double scale_tree_factor = getTreeScalingFactor(container, tree_);
    std::cout << maxState << std::endl;
    std::cout << countRange << std::endl;
    std::cout << scale_tree_factor*tree_->getTotalLength() << std::endl;
    tree_->scaleTree(scale_tree_factor);



    // Define substitution parameters
    std::map<int, std::vector<double>> paramMap = {
        {1, {-999}},
        {2, {3.0}},
        {3, {2.0, 0.1}},
        {4, {2.0, 0.1}},
        {5, {1.3}}
    };
    std::vector<int> rateChangeType = {8, 0, 1, 1, 0};
    int baseNum = -999;

    // Create substitution process
    auto mapOfNodeIds = getMapOfNodeIds(tree_);
    std::shared_ptr<DiscreteDistribution> gammaDist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree_);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree);
    

    // Create substitution model
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, paramMap, baseNum, countRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, false);
    // std::cout << chrModel->getBaseNumR()->getParameterValues()[0] << std::endl;
    subProcesses->addModel(chrModel, mapOfNodeIds[1]);
    SubstitutionProcess* nsubPro = subProcesses->clone();

    // Create likelihood object
    Context* context = new Context();
    bool weightedRootFreqs = true;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *container->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    

    std::vector<std::string> sequenceNames = container->getSequenceNames();
    std::cout << "Calculating likelihood" << std::endl;
    std::cout << newLik->getValue() << std::endl;
    auto substitutionModelParams = newLik->getSubstitutionModelParameters().getParameterNames();
    for (int i = 0; i < (int)(substitutionModelParams.size()); i++){
        std::cout << substitutionModelParams[i] << std::endl;
        std::cout << newLik->getParameters().getParameter(substitutionModelParams[i]).getValue() << std::endl;
    }
    std::cout << "Starting optimization" << std::endl;
    optimizeModelParametersOneDimension(newLik, 0.1, 2);
    std::cout << "Done" << std::endl;
    return 0;
}



VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    std::string line;
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::vector<std::string> familyIds;
    std::string columnName;

    // Skip the first two columns ("Organizem" and "Desc")
    std::getline(headerStream, columnName, '\t');  
    std::getline(headerStream, columnName, '\t');  

    while (std::getline(headerStream, columnName, '\t')) {
        familyIds.push_back(columnName);
    }

    // Create a sequence container
    auto* container = new VectorSiteContainer(alphabet);

    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string speciesName;
        std::vector<std::string> geneCounts;

        // Read species name
        std::getline(lineStream, speciesName, '\t');

        // Skip "Desc" column
        std::getline(lineStream, columnName, '\t');

        // Read gene counts
        while (std::getline(lineStream, columnName, '\t')) {
            geneCounts.push_back(columnName);
        }

        // Verify the number of gene counts matches the number of families
        if (geneCounts.size() != familyIds.size()) {
            delete container;
            throw std::runtime_error("Mismatch between family count and data columns in file.");
        }

        container->addSequence(BasicSequence(speciesName, geneCounts, alphabet), true);
    }
    file.close();
    return container;
}

std::vector<std::string> getSpeciesListFromFile(const std::string& filePath) {

    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    std::string line;

    if (!std::getline(file, line)) {
        throw std::runtime_error("File is empty or missing header: " + filePath);
    }

    // Parse the header
    std::istringstream headerStream(line);
    std::vector<std::string> speciesList;
    std::string columnName;

    // Skip the first two columns: "Desc" and "Family ID"
    std::getline(headerStream, columnName, '\t'); 
    std::getline(headerStream, columnName, '\t');

    // Read the remaining column headers (species names)
    while (std::getline(headerStream, columnName, '\t')) {
        speciesList.push_back(columnName);
    }

    file.close();
    return speciesList;
}

std::map<uint, vector<uint>> getMapOfNodeIds(PhyloTree* tree) {
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

std::tuple<int, int> getAlphabetLimit(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    std::string line;
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::string columnName;

    // Skip the first two columns ("Organizem" and "Desc")
    std::getline(headerStream, columnName, '\t');  
    std::getline(headerStream, columnName, '\t');  

    int max = 1;
    int min = 500;

    while (std::getline(file, line)) {
        std::string speciesName;
        std::istringstream lineStream(line);

        // Read species name
        std::getline(lineStream, speciesName, '\t');

        // Skip "Desc" column
        std::getline(lineStream, columnName, '\t');

        // Read gene counts
        while (std::getline(lineStream, columnName, '\t')) {
            int geneCount = std::stoi(columnName);
            if (geneCount > max) {
                max = geneCount;
            }
            if (geneCount < min) {
                min = geneCount;
            }
        }
    }
    file.close();
    return  std::make_tuple(min, max + 10);
}

int countUniqueStates(const Site site) {
    std::set<int> uniqueStates;
    // Iterate over all sequences (rows)
    for (size_t j = 0; j < site.size(); j++) {
        uniqueStates.insert(site[j]); 
    }
    return uniqueStates.size();
}

int getGeneCountRange(const VectorSiteContainer* container) {
    int min = 500;
    int max = 1;

    // Iterate over all sites (columns)
    for (size_t i = 0; i < container->getNumberOfSites(); i++) {
        const Site& site = container->getSite(i);
        for (size_t j = 0; j < site.size(); j++) {
            if (site[j] > max) {
                max = site[j];
            }
            if (site[j] < min) {
                min = site[j];
            }
        }
    }
    return max - min;
}

// Tol: 0.1, maxNumOfIterations: 2? 
void optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* likelihoodProcess, double tol, unsigned int maxNumOfIterations, bool mixed, unsigned curentIterNum)
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
    unsigned int numOfEvaluations = 0;
    int minDomain = 1; //TODO: hardcoded
    int maxDomain = 110; //TODO: hardcoded
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
    std::vector<int> rateChangeType = {8, 0, 1, 1, 0}; //TODO: hardcoded


    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model

    vector<string> parametersNames = likelihoodProcess->getSubstitutionModelParameters().getParameterNames();
    updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, parametersNames, 0, "");
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
            std::cout << "parameter value after optimization "+ std::to_string(likelihoodProcess->getParameters().getParameter(param.getName()).getValue()) << std::endl;
        }

        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();

    }
    delete optimizer;
}

void updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, std::vector<std::string> namesAllParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::string suffix){
    std::map<string, vector<std::pair<uint, int>>> sharedParamsNames;
    if (sharedParams){
        // LikelihoodUtils::createMapOfSharedParameterNames(*sharedParams, sharedParamsNames);
        std::cout << "Hi" << std::endl; 
    }
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
        if (sharedParams){
            if (sharedParamsNames.find(namesAllParams[i]) != sharedParamsNames.end()){
                for (size_t k = 0; k < sharedParamsNames[namesAllParams[i]].size(); k++ ){
                    uint model = sharedParamsNames[namesAllParams[i]][k].first;
                    int typeOfSharedParam = sharedParamsNames[namesAllParams[i]][k].second;
                    typeWithParamNames[typeOfSharedParam][model].push_back(namesAllParams[i]);

                }

            }

        }

    }

}

uint getModelFromParamName(string name){
    std::regex modelPattern ("_([\\d]+)");
    std::smatch sm;
    std::regex_search(name, sm, modelPattern);
    std::string modelSuffix = sm[sm.size()-1];
    uint modelId = static_cast<uint>(stoi(modelSuffix));
    return modelId;

}

void updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    
}

int getTypeOfParamFromParamName(string name){
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

double getTreeScalingFactor(const VectorSiteContainer* container, PhyloTree* tree) {
    int uniqueStateCount = 0;
    // Iterate over all sites (columns)
    for (size_t i = 0; i < container->getNumberOfSites(); i++) {
        const Site site = container->getSite(i);
        uniqueStateCount += countUniqueStates(site);
    }
    double avgUniqueStateCount = uniqueStateCount/container->getNumberOfSites();
    return avgUniqueStateCount/tree->getTotalLength();
}