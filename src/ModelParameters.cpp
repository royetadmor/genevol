#include "ModelParameters.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string> 
#include <tuple>
#include <stdexcept>



using namespace std;


ModelParameters::ModelParameters(BppApplication GenEvol)
{
    ModelParameters::treeFilePath_ = ApplicationTools::getAFilePath("_treePath", GenEvol.getParams(), true, true, "", true, "none", 1);
    ModelParameters::dataFilePath_ = ApplicationTools::getAFilePath("_dataPath", GenEvol.getParams(), true, true, "", true, "none", 1);
    ModelParameters::stateOverhead_ = ApplicationTools::getIntParameter("_stateOverhead", GenEvol.getParams(), STATE_OVERHEAD, "", true, -1);
    setAlphabetLimit(GenEvol);
    std::cout << ModelParameters::minState_ << ", " << ModelParameters::maxState_ << ", " << ModelParameters::stateOverhead_ << std::endl;
    ModelParameters::countRange_ = ApplicationTools::getIntParameter("_countRange", GenEvol.getParams(), ModelParameters::maxState_ - ModelParameters::minState_ + 1, "", true, 1);
    ModelParameters::alphabet_ = new IntegerAlphabet(ModelParameters::maxState_, ModelParameters::minState_);
    ModelParameters::container_ = readGeneFamilyFile(ModelParameters::dataFilePath_, ModelParameters::alphabet_);
    setBaseModelParameters(GenEvol);
    setRateFunctionTypes(GenEvol);

}

void ModelParameters::setAlphabetLimit(BppApplication GenEvol) {

    int minState = ApplicationTools::getIntParameter("_minState", GenEvol.getParams(), -1, "", true, 1);
    int maxState = ApplicationTools::getIntParameter("_maxState", GenEvol.getParams(), -1, "", true, 1);

    // Both are given by the user
    if (minState != -1 && maxState != -1) {
        ModelParameters::minState_ = minState;
        ModelParameters::maxState_ = maxState;
        return;
    }
    std::ifstream file(ModelParameters::dataFilePath_);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + ModelParameters::dataFilePath_);
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
    if (minState == -1) {
        ModelParameters::minState_ = min;    
    }
    if (maxState == -1) {
        ModelParameters::maxState_ = max + ModelParameters::stateOverhead_;    
    }
}

void ModelParameters::setBaseModelParameters(BppApplication GenEvol) {
    ModelParameters::paramMap_ = {
        {1, ApplicationTools::getVectorParameter<double>("_baseNumR", GenEvol.getParams(), ',', "-999", "", true, -1)},
        {2, ApplicationTools::getVectorParameter<double>("_dupl", GenEvol.getParams(), ',', "1", "", true, -1)},
        {3, ApplicationTools::getVectorParameter<double>("_loss", GenEvol.getParams(), ',', "1,0.1", "", true, -1)},
        {4, ApplicationTools::getVectorParameter<double>("_gain", GenEvol.getParams(), ',', "1,0.1", "", true, -1)},
        {5, ApplicationTools::getVectorParameter<double>("_demi", GenEvol.getParams(), ',', "-999", "", true, -1)}
    };
}

void ModelParameters::setRateFunctionTypes(BppApplication GenEvol) {
    const int gainFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_gainFunc", GenEvol.getParams(), "LINEAR", "", true, -1));
    const int lossFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_lossFunc", GenEvol.getParams(), "LINEAR", "", true, -1));
    const int duplFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_duplFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int demiDuplFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_demiDuplFunc", GenEvol.getParams(), "IGNORE", "", true, -1));
    const int baseNumRFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_baseNumRFunc", GenEvol.getParams(), "IGNORE", "", true, -1));
    ModelParameters::rateChangeType_.push_back(baseNumRFunc);
    ModelParameters::rateChangeType_.push_back(duplFunc);
    ModelParameters::rateChangeType_.push_back(lossFunc);
    ModelParameters::rateChangeType_.push_back(gainFunc);
    ModelParameters::rateChangeType_.push_back(demiDuplFunc);
}

VectorSiteContainer* ModelParameters::readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet) {
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