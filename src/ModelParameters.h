#ifndef GENEVOL_MODELPARAMETERS_H
#define GENEVOL_MODELPARAMETERS_H

#define STATE_OVERHEAD 10

#include <string> 
#include <tuple>

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>


#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>



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
    IntegerAlphabet* alphabet_;
    VectorSiteContainer* container_;
    std::map<int, std::vector<double>> paramMap_;
    std::vector<int> rateChangeType_;
    const std::map<std::string, int> func_string_to_enum = {
        {"CONST", 0},
        {"LINEAR", 1},
        {"IGNORE", 8} // Can add more according to enum at ChromosomeSubModel.h
    };
    
public:
    ModelParameters(BppApplication GenEvol);
    ~ModelParameters(){};

private:
    void setAlphabetLimit(BppApplication GenEvol);
    void setBaseModelParameters(BppApplication GenEvol);
    VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
    void setRateFunctionTypes(BppApplication GenEvol);
};


#endif // GENEVOL_MODELPARAMETERS_H