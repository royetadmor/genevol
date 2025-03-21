#include <iostream>
#include <string> 
#include <set>

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
double countUniqueStates(const VectorSiteContainer* container);
VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
std::vector<std::string> getSpeciesListFromFile(const std::string& filePath);
std::map<uint, vector<uint>> getMapOfNodeIds(PhyloTree* tree);


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
    double uniqueStatesCount = countUniqueStates(container);
    int countRange = getGeneCountRange(container);
    std::cout << maxState << std::endl;
    std::cout << uniqueStatesCount << std::endl;
    std::cout << countRange << std::endl;
    double treeLength = tree_->getTotalLength();

    // TODO: hardcoded
    auto scale_tree_factor = 2;//uniqueStatesCount/treeLength; // rescale by average amount of unique state by position
    countRange = 50;
    tree_->scaleTree(scale_tree_factor);



    // Define substitution parameters
    std::map<int, std::vector<double>> paramMap = {
        {1, {2.0}},
        {2, {3.0}},
        {3, {2.0, 0.1}},
        {4, {2.0, 0.1}},
        {5, {1.3}}
    };
    std::vector<int> rateChangeType = {0, 0, 1, 1, 0};
    int baseNum = 7;

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
    std::cout << container->getSequence(sequenceNames[1]).toString() << std::endl;
    std::cout << "Calculating likelihood" << std::endl;
    std::cout << newLik->getValue() << std::endl;
    auto substitutionModelParams = newLik->getSubstitutionModelParameters().getParameterNames();
    for (int i = 0; i < (int)(substitutionModelParams.size()); i++){
        std::cout << substitutionModelParams[i] << std::endl;
       std::cout << newLik->getParameters().getParameter(substitutionModelParams[i]).getValue() << std::endl;
    }
    // std::cout << substitutionModelParams.size() << std::endl;
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

double countUniqueStates(const VectorSiteContainer* container) {
    std::set<int> uniqueStates;
    // Iterate over all sites (columns)
    for (size_t i = 0; i < container->getNumberOfSites(); i++) {
        const Site& site = container->getSite(i);

        // Iterate over all sequences (rows)
        for (size_t j = 0; j < site.size(); j++) {
            uniqueStates.insert(site[j]); 
        }
    }
    return static_cast<double>(uniqueStates.size());
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
