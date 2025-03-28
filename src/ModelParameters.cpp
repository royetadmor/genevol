#include "ModelParameters.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string> 
#include <tuple>



using namespace std;


ModelParameters::ModelParameters(const string& treePath, const string& dataPath)
{
    ModelParameters::treeFilePath_ = treePath;
    ModelParameters::dataFilePath_ = dataPath;
    setAlphabetLimit();
    ModelParameters::countRange_ = ModelParameters::maxState_ - 10 - ModelParameters::minState_;

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
    ModelParameters::minState_ = min;
    ModelParameters::maxState_ = max + 10;
}