

#ifndef GENECOUNT_DEPENDENCY_FUNC_H
#define GENECOUNT_DEPENDENCY_FUNC_H

// From bpp-core
#include <Bpp/Numeric/Function/Functions.h>


#define lowerBoundOfRateParam 0.0
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0
#define lowerBoundOfExpParam -3.0
#define upperBoundExpParam 4.6
#define revSigmoidExpRateParam 2
#define logNormalDomainFactor 4


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

    class ConstantDependencyFunction :
        public virtual GeneCountDependencyFunction
    {
    public:
        ConstantDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~ConstantDependencyFunction(){}

        FunctionType getName() const{return FunctionType::CONSTANT;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 1;}

    };

    class LinearDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:

        LinearDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~LinearDependencyFunction(){}

        FunctionType getName() const{return FunctionType::LINEAR;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 2;}
        void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum);
        void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound);
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
    };

    class LinearBDDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        LinearBDDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~LinearBDDependencyFunction(){}

        FunctionType getName() const{return FunctionType::LINEAR_BD;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 1;}

    };
    
    class ExponentailDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        ExponentailDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~ExponentailDependencyFunction(){}

        FunctionType getName() const {return FunctionType::EXP;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 2;}
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    };

    class PolynomialDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        PolynomialDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~PolynomialDependencyFunction(){}

        FunctionType getName() const{return FunctionType::POLYNOMIAL;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 3;}
        void setDomainsIfNeeded(int minChrNum, int maxChrNum)
        {
            domainMin_ = minChrNum;
            domainMax_ = maxChrNum;
        }
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    };

    class LognormalDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        LognormalDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~LognormalDependencyFunction(){}

        FunctionType getName() const {return FunctionType::LOGNORMAL;}
        void setDomainsIfNeeded(int minChrNum, int maxChrNum)
        {
            domainMin_ = minChrNum;
            domainMax_ = maxChrNum;
        }

        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 3;}
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    };

    class RevSigmoidDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        RevSigmoidDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~RevSigmoidDependencyFunction(){}

        FunctionType getName() const {return FunctionType::REVERSE_SIGMOID;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 3;}
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);
        void setDomainsIfNeeded(int minChrNum, int maxChrNum) {
            domainMin_ = minChrNum;
            domainMax_ = maxChrNum;
        }

    };

    class LogitnormalDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        LogitnormalDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~LogitnormalDependencyFunction(){}

        FunctionType getName() const {return FunctionType::LOGITNORMAL;}
        void setDomainsIfNeeded(int minChrNum, int maxChrNum)
        {
            domainMin_ = minChrNum;
            domainMax_ = maxChrNum;
        }

        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 3;}
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    };

}


#endif // GENECOUNT_DEPENDENCY_FUNC_H