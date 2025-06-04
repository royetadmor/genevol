#include "ModelParameters.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string> 
#include <tuple>
#include <stdexcept>



using namespace std;


ModelParameters::ModelParameters()
{
    ModelParameters::treeFilePath_ = ModelParameters::getEnvVar("TREE_PATH");
    ModelParameters::dataFilePath_ = ModelParameters::getEnvVar("DATA_PATH");
    setAlphabetLimit();
    ModelParameters::countRange_ = ModelParameters::getIntVar("COUNT_RANGE", ModelParameters::maxState_ - 10 - ModelParameters::minState_);
    ModelParameters::alphabet_ = new IntegerAlphabet(ModelParameters::maxState_, ModelParameters::minState_);
    ModelParameters::container_ = readGeneFamilyFile(ModelParameters::dataFilePath_, ModelParameters::alphabet_);
    setBaseModelParameters();
    setRateFunctionTypes();

}

void ModelParameters::setAlphabetLimit() {
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
    ModelParameters::minState_ = ModelParameters::getIntVar("MIN_STATE", min);
    ModelParameters::maxState_ = ModelParameters::getIntVar("MAX_STATE", max + 10);
}

void ModelParameters::setBaseModelParameters() {
    ModelParameters::paramMap_ = {
        {1, ModelParameters::getVectorVar("BASE_NUM_R", -999)},
        {2, ModelParameters::getVectorVar("DUPL", 0)},
        {3, ModelParameters::getVectorVar("GAIN", 0)},
        {4, ModelParameters::getVectorVar("LOSS", 0)},
        {5, ModelParameters::getVectorVar("DEMI", 0)}
    };
}

void ModelParameters::setRateFunctionTypes() {
    const int gainFunc = ModelParameters::func_string_to_enum.at(ModelParameters::getEnvVar("GAIN_FUNC"));
    const int lossFunc = ModelParameters::func_string_to_enum.at(ModelParameters::getEnvVar("LOSS_FUNC"));
    const int duplFunc = ModelParameters::func_string_to_enum.at(ModelParameters::getEnvVar("DUPL_FUNC"));
    const int demiDuplFunc = ModelParameters::func_string_to_enum.at(ModelParameters::getEnvVar("DEMI_DUPL_FUNC"));
    const int baseNumRFunc = ModelParameters::func_string_to_enum.at(ModelParameters::getEnvVar("BASE_NUM_R_FUNC"));
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

std::string ModelParameters::getEnvVar(const std::string& key) {
    const char* val = std::getenv(key.c_str());
    if (val == nullptr) {
        const std::string err = "Failed to find parameter " + key +". Exiting.";
        throw Exception(err);
    }
    return std::string(val);
}

int ModelParameters::getIntVar(const std::string& key, const int defaultVal) {
    int val;
    try {
        val = std::stoi(ModelParameters::getEnvVar(key));
    } catch (bpp::Exception& e) {
        val = defaultVal;
    }
    return val;
}

std::vector<double> ModelParameters::getVectorVar(const std::string& key, const double defaultVal) {
    std::vector<double> val;
    try {
        const std::string envVarValue = ModelParameters::getEnvVar(key);
        std::stringstream ss(envVarValue);
        std::string token;
        while (std::getline(ss, token, ',')) {
            try {
                val.push_back(std::stod(token));
            } catch (const std::invalid_argument& e) {
                throw std::runtime_error("Invalid double value: '" + token + "'");
            }
        }

    } catch (bpp::Exception& e) {
        std::cout << "Caught an exception while parsing " + key << std::endl;
        std::cout << e.what() << std::endl;
        val.push_back(defaultVal);
    }
    return val;
}