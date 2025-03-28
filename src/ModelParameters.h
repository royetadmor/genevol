#include <string> 
#include <tuple>


using namespace std;

class ModelParameters
{

public: // Variables
    string treeFilePath_;
    string dataFilePath_;
    int minState_;
    int maxState_;
    int countRange_;
    
public:
    ModelParameters(const string& treePath, const string& dataPath);
    ~ModelParameters(){};

    void setAlphabetLimit();
};
