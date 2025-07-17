#ifndef GENEVOL_MODELPARAMETERS_H
#define GENEVOL_MODELPARAMETERS_H

#define STATE_OVERHEAD 10

#include <string> 
#include <tuple>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>


#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include "ChromosomeSubstitutionModel.h"

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
    IntegerAlphabet* alphabet_;
    VectorSiteContainer* container_;
    std::map<int, std::vector<double>> paramMap_;
    std::vector<int> rateChangeType_;

    // Translates rate functions from string to int
    const std::map<std::string, int> func_string_to_enum = {
        {"CONST", ChromosomeNumberDependencyFunction::FunctionType::CONSTANT},
        {"LINEAR", ChromosomeNumberDependencyFunction::FunctionType::LINEAR},
        {"IGNORE", ChromosomeNumberDependencyFunction::FunctionType::IGNORE} // Can add more according to enum at ChromosomeSubModel.h
    };

    // Expected number of parameters for each rate function
    const std::map<int, int> expectedNumOfParams = {
        {ChromosomeNumberDependencyFunction::FunctionType::CONSTANT, 1},
        {ChromosomeNumberDependencyFunction::FunctionType::LINEAR, 2},
        {ChromosomeNumberDependencyFunction::FunctionType::IGNORE, 1}
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
};


#endif // GENEVOL_MODELPARAMETERS_H