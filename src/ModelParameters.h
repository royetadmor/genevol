#ifndef GENEVOL_MODELPARAMETERS_H
#define GENEVOL_MODELPARAMETERS_H

#include <string> 
#include <tuple>


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
    ModelParameters();
    ~ModelParameters(){};

private:
    void setAlphabetLimit();
    void setBaseModelParameters();
    VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
    std::string getEnvVar(const std::string& key);
    int getIntVar(const std::string& key, const int defaultVal);
    std::vector<double> getVectorVar(const std::string& key, const double defaultVal);
    void setRateFunctionTypes();
};


#endif // GENEVOL_MODELPARAMETERS_H