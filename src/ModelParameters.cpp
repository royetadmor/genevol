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
    // Load general params
    ModelParameters::treeFilePath_ = ApplicationTools::getAFilePath("_treePath", GenEvol.getParams(), true, true, "", true, "none", 1);
    ModelParameters::dataFilePath_ = ApplicationTools::getAFilePath("_dataPath", GenEvol.getParams(), true, true, "", true, "none", 1);
    ModelParameters::stateOverhead_ = ApplicationTools::getIntParameter("_stateOverhead", GenEvol.getParams(), STATE_OVERHEAD, "", true, -1);
    setAlphabetLimit(GenEvol);
    ModelParameters::alphabet_ = std::make_shared<bpp::GeneCountAlphabet>(ModelParameters::maxState_, ModelParameters::minState_);
    ModelParameters::container_ = readGeneFamilyFile(ModelParameters::dataFilePath_, ModelParameters::alphabet_);
    ModelParameters::rDist_ = PhylogeneticsApplicationTools::getRateDistribution(GenEvol.getParams(), "", true, true);

    // Load customization params
    ModelParameters::countRange_ = ApplicationTools::getIntParameter("_countRange", GenEvol.getParams(), ModelParameters::maxState_ - ModelParameters::minState_ + 1, "", true, 1);
    ModelParameters::branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", GenEvol.getParams(), -999.0);

    // Load MM params
    ModelParameters::alphaGain_ = ApplicationTools::getDoubleParameter("_alphaGain", GenEvol.getParams(), 1.0);
    ModelParameters::betaGain_ = ApplicationTools::getDoubleParameter("_betaGain", GenEvol.getParams(), 1.0);
    ModelParameters::alphaLoss_ = ApplicationTools::getDoubleParameter("_alphaLoss", GenEvol.getParams(), 1.0);
    ModelParameters::betaLoss_ = ApplicationTools::getDoubleParameter("_betaLoss", GenEvol.getParams(), 1.0);
    ModelParameters::categories_ = ApplicationTools::getIntParameter("_numCategories", GenEvol.getParams(), 4, "", true, -1);
    ModelParameters::mixtureInnovation_ = ApplicationTools::getDoubleParameter("_mixtureInnovation", GenEvol.getParams(), 1.0);
    ModelParameters::mixtureElimination_ = ApplicationTools::getDoubleParameter("_mixtureElimination", GenEvol.getParams(), 1.0);

    // Load parameters and rate functions
    setBaseModelParameters(GenEvol);
    setRateFunctionTypes(GenEvol);
    validateRateFunctionParameters();

    // Set fixed params
    ModelParameters::fixedParams_ = ApplicationTools::getVectorParameter<string>("_fixedParams", GenEvol.getParams(), ',', "", "", true, 1);
    ModelParameters::mixtureFixedParams_ = ApplicationTools::getVectorParameter<string>("_mixtureFixedParams", GenEvol.getParams(), ',', "", "", true, 1);

    // Set constrainted params
    std::vector<string> constraintedParams = ApplicationTools::getVectorParameter<string>("_constraintedParams", GenEvol.getParams(), ';', "", "", true, 1);
    setConstraintedParams(GenEvol, constraintedParams, ModelParameters::constraintedParams_);
    std::vector<string> mixtureConstraintedParams = ApplicationTools::getVectorParameter<string>("_mixtureConstraintedParams", GenEvol.getParams(), ';', "", "", true, 1);
    setConstraintedParams(GenEvol, mixtureConstraintedParams, ModelParameters::mixtureConstraintedParams_);
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
        {0, ApplicationTools::getVectorParameter<double>("_loss", GenEvol.getParams(), ',', "1", "", true, -1)},
        {1, ApplicationTools::getVectorParameter<double>("_gain", GenEvol.getParams(), ',', "1", "", true, -1)},
        {2, ApplicationTools::getVectorParameter<double>("_innovation", GenEvol.getParams(), ',', "1", "", true, -1)},
        {3, ApplicationTools::getVectorParameter<double>("_elimination", GenEvol.getParams(), ',', "1", "", true, -1)}
    };
}

void ModelParameters::setRateFunctionTypes(BppApplication GenEvol) {
    const int gainFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_gainFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int lossFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_lossFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int innovationFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_innovationFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int eliminationFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_eliminationFunc", GenEvol.getParams(), "CONST", "", true, -1));

    const int mixtureGainFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_mixtureGainFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int mixtureLossFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_mixtureLossFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int mixtureInnovationFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_mixtureInnovationFunc", GenEvol.getParams(), "CONST", "", true, -1));
    const int mixtureEliminationFunc = ModelParameters::func_string_to_enum.at(ApplicationTools::getStringParameter("_mixtureEliminationFunc", GenEvol.getParams(), "CONST", "", true, -1));

    ModelParameters::rateChangeType_.push_back(lossFunc);
    ModelParameters::rateChangeType_.push_back(gainFunc);
    ModelParameters::rateChangeType_.push_back(innovationFunc);
    ModelParameters::rateChangeType_.push_back(eliminationFunc);

    ModelParameters::mixtureRateChangeType_.push_back(mixtureLossFunc);
    ModelParameters::mixtureRateChangeType_.push_back(mixtureGainFunc);
    ModelParameters::mixtureRateChangeType_.push_back(mixtureInnovationFunc);
    ModelParameters::mixtureRateChangeType_.push_back(mixtureEliminationFunc);
}

void ModelParameters::validateRateFunctionParameters() {
    // For each rate category (in order)
    for (size_t i = 0; i < ModelParameters::rateChangeType_.size(); ++i) {
        int funcType = rateChangeType_[i];
        int paramCategoryId = static_cast<int>(i);

        auto it = expectedNumOfParams.find(funcType);
        if (it == expectedNumOfParams.end()) {
            throw std::runtime_error("Unknown function type enum: " + std::to_string(funcType));
        }

        size_t requiredParams = it->second;
        size_t actualParams = ModelParameters::paramMap_.at(paramCategoryId).size();

        if (actualParams != requiredParams) {
            throw std::runtime_error(
                "Parameter validation failed for category id " + std::to_string(paramCategoryId) +
                ": function type requires " + std::to_string(requiredParams) +
                " parameters, but got " + std::to_string(actualParams)
            );
        }
    }
}

std::shared_ptr<VectorSiteContainer> ModelParameters::readGeneFamilyFile(const std::string& filePath, std::shared_ptr<const bpp::Alphabet> alphabet) {
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
    auto container = std::make_shared<VectorSiteContainer>(alphabet);

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
            throw std::runtime_error("Mismatch between family count and data columns in file.");
        }

        std::unique_ptr<bpp::Sequence> seq = std::make_unique<bpp::Sequence>(speciesName, geneCounts, alphabet);
        container->addSequence(speciesName, seq);
    }
    file.close();
    return container;
}

void ModelParameters::setConstraintedParams(BppApplication GenEvol, std::vector<string> inputParams, std::map<string, vector<string>>& outputParams) {
    for (const string& element: inputParams) {
        std::vector<std::string> paramList;
        std::stringstream ss(element); 
        std::string segment;

        // split string by comma
        while (std::getline(ss, segment, ',')) {
            paramList.push_back(segment);
        }
        if (paramList.size() != 2) {
            throw std::runtime_error("Each constraint should contain exactly 2 parameters");
        }
        outputParams[paramList[0]].push_back(paramList[1]);
    }
}