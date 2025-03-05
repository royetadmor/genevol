#include <iostream>
#include <string> 


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



VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
std::vector<std::string> getSpeciesListFromFile(const std::string& filePath);
std::map<uint, vector<uint>> getMapOfNodeIds(PhyloTree* tree);


int main(int args, char **argv) {
    // Get tree path and count path from args
    IntegerAlphabet* alphabet = new IntegerAlphabet(500); // Max gene family members, wild guess
    std::string treePath = "/workspaces/genevol/tree.newick";
    std::string dataPath = "/workspaces/genevol/data.fasta";

    // Get tree
    Newick reader;
    PhyloTree* tree_ = reader.readPhyloTree(treePath);

    VectorSiteContainer* container = readGeneFamilyFile(dataPath, alphabet);
    std::vector<std::string> species = getSpeciesListFromFile(dataPath);

    // Define substitution parameters
    std::map<int, std::vector<double>> paramMap = {
        {1, {2.0}},
        {2, {3.0}},
        {3, {2.0}},
        {4, {2.0}},
        {5, {1.3}}
    };
    std::vector<int> rateChangeType = {0, 0, 0, 0, 0};
    std::map<uint, std::vector<int>> fixedParams = {
        {1, {1, 1, 1, 1, 1}} 
    };
    int baseNum = 70;

    // Create substitution process
    auto mapOfNodeIds = getMapOfNodeIds(tree_);
    std::shared_ptr<DiscreteDistribution> gammaDist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    auto parTree = std::make_shared<ParametrizablePhyloTree>(tree_);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProcesses = std::make_shared<NonHomogeneousSubstitutionProcess>(gammaDist, parTree);
    

    // Create substitution model
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, paramMap, baseNum, 500, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChangeType, false, false);
    subProcesses->addModel(chrModel, mapOfNodeIds[1]);
    SubstitutionProcess* nsubPro = subProcesses->clone();

    // Create likelihood object
    Context* context = new Context();
    bool weightedRootFreqs = true;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *container->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    

    std::vector<std::string> sequenceNames = container->getSequenceNames();
    std::cout << newLik->getValue() << std::endl;
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

        container->addSequence(BasicSequence(speciesName, geneCounts[0], alphabet), true);
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
