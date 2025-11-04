#ifndef GENEVOL_MODELPARAMETERS_H
#define GENEVOL_MODELPARAMETERS_H

#define STATE_OVERHEAD 10

#include <string> 
#include <tuple>
#include <sstream>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>


#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

#include "GeneCountDependencyFunction.h"

using namespace std;
using namespace bpp;

class ModelParameters
{
public: // Variables
    string treeFilePath_;
    string dataFilePath_;
    int minState_;
    int maxState_;
    int stateOverhead_;
    int countRange_;
    int categories_;
    double branchMul_;
    double alphaGain_;
    double betaGain_;
    double alphaLoss_;
    double betaLoss_;
    double mixtureInnovation_;
    double mixtureElimination_;
    IntegerAlphabet* alphabet_;
    VectorSiteContainer* container_;
    std::map<int, std::vector<double>> paramMap_;
    std::vector<int> mixtureRateChangeType_;
    std::vector<int> rateChangeType_;
    std::vector<string> fixedParams_;
    std::vector<string> mixtureFixedParams_;
    std::map<string, vector<string>> constraintedParams_;
    std::map<string, vector<string>> mixtureConstraintedParams_;
    DiscreteDistribution* rDist_;

    // Translates rate functions from string to int
    const std::map<std::string, int> func_string_to_enum = {
        {"CONST", GeneCountDependencyFunction::FunctionType::CONSTANT},
        {"LINEAR", GeneCountDependencyFunction::FunctionType::LINEAR},
        {"EXP", GeneCountDependencyFunction::FunctionType::EXP},
        {"IGNORE", GeneCountDependencyFunction::FunctionType::IGNORE} // Can add more according to enum at GeneCountSubModel.h
    };

    // Expected number of parameters for each rate function
    const std::map<int, int> expectedNumOfParams = {
        {GeneCountDependencyFunction::FunctionType::CONSTANT, 1},
        {GeneCountDependencyFunction::FunctionType::LINEAR, 2},
        {GeneCountDependencyFunction::FunctionType::EXP, 2},
        {GeneCountDependencyFunction::FunctionType::IGNORE, 1}
    };
    
public:
    ModelParameters(BppApplication GenEvol);
    ~ModelParameters(){};

private:
    void setAlphabetLimit(BppApplication GenEvol);
    void setBaseModelParameters(BppApplication GenEvol);
    VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
    void setRateFunctionTypes(BppApplication GenEvol);
    void validateRateFunctionParameters();
    void setConstraintedParams(BppApplication GenEvol, std::vector<string> inputParams, std::map<string, vector<string>>& outputParams);
};


#endif // GENEVOL_MODELPARAMETERS_H