

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
    };

    class NLinearBDDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NLinearBDDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NLinearBDDependencyFunction(){}

        FunctionType getName() const{return FunctionType::LINEAR_BD;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 1;}

    };
    
    class NExponentailDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NExponentailDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NExponentailDependencyFunction(){}

        FunctionType getName() const {return FunctionType::EXP;}
        double getRate(std::vector<Parameter*> params, size_t state) const;
        size_t getNumOfParameters() const{return 2;}
        void getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber);
        void getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    };

    class NPolynomialDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NPolynomialDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NPolynomialDependencyFunction(){}

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

    class NLognormalDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NLognormalDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NLognormalDependencyFunction(){}

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

    class NRevSigmoidDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NRevSigmoidDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NRevSigmoidDependencyFunction(){}

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

    class NLogitnormalDependencyFunction:
        public virtual GeneCountDependencyFunction
    {
    public:
        NLogitnormalDependencyFunction():GeneCountDependencyFunction(){}
        virtual ~NLogitnormalDependencyFunction(){}

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