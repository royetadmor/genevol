#ifndef GENEVOL_MODELPARAMETERS_H
#define GENEVOL_MODELPARAMETERS_H

#define STATE_OVERHEAD 10
#define STATE_SPACE_UPPER_BOUND 200

#include <string> 
#include <tuple>
#include <sstream>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>


#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

#include "GeneCountDependencyFunction.h"
#include "GeneCountAlphabet.h"

using namespace std;
using namespace bpp;

class ModelParameters
{
public: // Variables
    string treeFilePath_;
    string dataFilePath_;
    int minState_;
    int maxState_;
    bool allowCapState_;
    int stateOverhead_;
    int categories_;
    double branchMul_;
    double alphaGain_;
    double betaGain_;
    double alphaLoss_;
    double betaLoss_;
    double mixtureInnovation_;
    double mixtureElimination_;
    std::shared_ptr<GeneCountAlphabet> alphabet_;
    std::shared_ptr<VectorSiteContainer> container_;
    std::map<int, std::vector<double>> paramMap_;
    std::vector<int> mixtureRateChangeType_;
    std::vector<int> rateChangeType_;
    std::vector<string> fixedParams_;
    std::vector<string> mixtureFixedParams_;
    std::map<string, string> constraintedParams_;
    std::map<string, string> mixtureConstraintedParams_;
    std::shared_ptr<DiscreteDistributionInterface> rDist_;

    // Translates rate functions from string to int
    const std::map<std::string, int> func_string_to_enum = {
        {"CONST", GeneCountDependencyFunction::FunctionType::CONSTANT},
        {"LINEAR", GeneCountDependencyFunction::FunctionType::LINEAR},
        {"LINEAR_BD", GeneCountDependencyFunction::FunctionType::LINEAR_BD},
        {"EXP", GeneCountDependencyFunction::FunctionType::EXP},
        {"POLYNOMIAL", GeneCountDependencyFunction::FunctionType::POLYNOMIAL},
        {"LOGNORMAL", GeneCountDependencyFunction::FunctionType::LOGNORMAL},
        {"REV_SIGMOID", GeneCountDependencyFunction::FunctionType::REVERSE_SIGMOID},
        // {"LOGITNORMAL", GeneCountDependencyFunction::FunctionType::LOGITNORMAL}, TODO: this currently doesn't work
        {"IGNORE", GeneCountDependencyFunction::FunctionType::IGNORE}
    };

    // Expected number of parameters for each rate function
    const std::map<int, int> expectedNumOfParams = {
        {GeneCountDependencyFunction::FunctionType::CONSTANT, 1},
        {GeneCountDependencyFunction::FunctionType::LINEAR, 2},
        {GeneCountDependencyFunction::FunctionType::LINEAR_BD, 1},
        {GeneCountDependencyFunction::FunctionType::EXP, 2},
        {GeneCountDependencyFunction::FunctionType::POLYNOMIAL, 3},
        {GeneCountDependencyFunction::FunctionType::LOGNORMAL, 3},
        {GeneCountDependencyFunction::FunctionType::REVERSE_SIGMOID, 3},
        // {GeneCountDependencyFunction::FunctionType::LOGITNORMAL, 3}, TODO: this currently doesn't work
        {GeneCountDependencyFunction::FunctionType::IGNORE, 1}
    };
    
public:
    ModelParameters(BppApplication GenEvol);
    ~ModelParameters(){};

private:
    void setAlphabetLimit(BppApplication GenEvol);
    void setBaseModelParameters(BppApplication GenEvol);
    std::shared_ptr<VectorSiteContainer> readGeneFamilyFile(const std::string& filePath, std::shared_ptr<const bpp::Alphabet> alphabet);
    void setRateFunctionTypes(BppApplication GenEvol);
    void validateRateFunctionParameters();
    void setConstraintedParams(BppApplication GenEvol, std::vector<string> inputParams, std::map<string, string>& outputParams);
    std::string capState(std::string geneCount);
};


#endif // GENEVOL_MODELPARAMETERS_H