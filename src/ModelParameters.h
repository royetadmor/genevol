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
    
public:
    ModelParameters(const string& treePath, const string& dataPath);
    ~ModelParameters(){};

private:
    void setAlphabetLimit();
    VectorSiteContainer* readGeneFamilyFile(const std::string& filePath, IntegerAlphabet* alphabet);
};
