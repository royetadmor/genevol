

#ifndef GENECOUNT_DEPENDENCY_FUNC_H
#define GENECOUNT_DEPENDENCY_FUNC_H

// From bpp-core
#include <Bpp/Numeric/Function/Functions.h>


#define lowerBoundOfRateParam 0.0
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0


using namespace std;

namespace bpp
{
    class GeneCountDependencyFunction
    {
    protected:
        // the min and max chromosome counts
        int domainMin_;
        int domainMax_;

    public:
        enum FunctionType {CONSTANT, LINEAR, LINEAR_BD, EXP, POLYNOMIAL, LOGNORMAL, REVERSE_SIGMOID, LOGITNORMAL, IGNORE};
        GeneCountDependencyFunction():domainMin_(0), domainMax_(0){}
        virtual ~GeneCountDependencyFunction(){}

        virtual FunctionType getName() const = 0;
        virtual double getRate(std::vector<Parameter*> params, size_t state) const = 0;
        virtual size_t getNumOfParameters() const = 0;
        virtual void setDomainsIfNeeded(int minChrNum, int maxChrNum){}

        virtual void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum){
            std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
            *lowerBound = interval->getLowerBound();
            *upperBound = interval->getUpperBound();
        }

        virtual void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){return;};
        virtual void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber)
        {
            *lowerBound = lowerBoundOfRateParam;
            *upperBound = upperBoundOfRateParam;
        }
        virtual void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
            *lowerBound = lowerBoundOfRateParam;
            *upperBound = upperBoundOfRateParam;
        }
    };

    class NConstantDependencyFunction : //TODO: change after migration
        public virtual GeneCountDependencyFunction
    {
    public:
        NConstantDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NConstantDependencyFunction(){}

        FunctionType getName() const{return FunctionType::CONSTANT;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 1;}

    };

    class NLinearDependencyFunction: //TODO: change after migration
        public virtual GeneCountDependencyFunction
    {
    public:

        NLinearDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NLinearDependencyFunction(){}

        FunctionType getName() const{return FunctionType::LINEAR;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 2;}
        void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum);
        void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
        double getParsimonyBound(std::vector<double> params, double parsimonyBound, size_t index, int minChrNum, int maxChrNum);

    };
}


#endif // GENECOUNT_DEPENDENCY_FUNC_H