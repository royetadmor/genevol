#include <iostream>
#include <string> 


#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Io/Newick.h>


using namespace bpp;
using namespace std;



VectorSequenceContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
int validateTreeSpecies(const PhyloTree& tree, const std::vector<std::string>& speciesList);
std::vector<std::string> getSpeciesListFromFile(const std::string& filePath);


int main(int args, char **argv) {
    // Get tree path and count path from args
    IntegerAlphabet* alphabet = new IntegerAlphabet(500); // Max gene family members, wild guess
    std::string treePath = "/workspaces/genevol/tree.newick";
    std::string dataPath = "/workspaces/genevol/data.fasta";

    // Get tree
    Newick reader;
    PhyloTree* tree_ = reader.readPhyloTree(treePath);

    VectorSequenceContainer* container = readGeneFamilyFile(dataPath, alphabet);
    std::vector<std::string> species = getSpeciesListFromFile(dataPath);
    
    std::cout << container->getSequence(2).toString() << std::endl;
    std::cout << "Done" << std::endl;
    return 0;
}

VectorSequenceContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet) {

    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filePath);
    }

    // Read the header to extract species names
    std::string line;
    std::getline(file, line);
    std::istringstream headerStream(line);
    std::vector<std::string> speciesNames;
    std::string columnName;

    // Skip the first two columns ("Desc" and "Family ID")
    std::getline(headerStream, columnName, '\t');
    std::getline(headerStream, columnName, '\t');

    while (std::getline(headerStream, columnName, '\t')) {
        speciesNames.push_back(columnName);
    }

    auto* container = new VectorSequenceContainer(alphabet);
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string familyId;
        std::vector<std::string> geneCounts;

        // Skip the "Desc" column
        std::getline(lineStream, columnName, '\t');

        // Read the Family ID
        std::getline(lineStream, familyId, '\t');

        // Read the gene counts
        while (std::getline(lineStream, columnName, '\t')) {
            geneCounts.push_back(columnName);
        }

        // Verify the number of gene counts matches the number of species
        if (geneCounts.size() != speciesNames.size()) {
            delete container;
            throw std::runtime_error("Mismatch between species count and data columns in file.");
        }

        auto* seq = new BasicSequence(familyId, geneCounts, alphabet);
        container->addSequence(*seq, true);
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
