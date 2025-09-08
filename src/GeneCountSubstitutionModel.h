#ifndef GENECOUNT_SUBSTITUTION_MODEL_H
#define GENECOUNT_SUBSTITUTION_MODEL_H

// from bpp-phyl
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Phyl/Model/AbstractSubstitutionModel.h>

// from bpp-core
#include <Bpp/Numeric/Matrix/MatrixTools.h>

// from bpp-seq
#include <Bpp/Seq/Alphabet/IntegerAlphabet.h>

// local
#include "GeneCountDependencyFunction.h"
#include "ModelParameters.h"

#define IgnoreParam -999
#define NUM_MODEL_PARAMETERS 4
#define TAYLOR_ITERATION_EPSILON 0.0001

using namespace std;

namespace bpp
{
  class GeneCountSubstitutionModel;

  class NcompositeParameter{ //TODO: change naming after migration

  private:
    // a vector of ChromosomeSubstitutionModel parameters that correspond to a given parameter type
    std::vector<Parameter*> params_; 

    // A function which is used to calculate the bounds and the rate of the composite parameter
    GeneCountDependencyFunction* func_; 

    // The name of the rate parameter. For exampe, "gain" for a paremter that contains gain0, gain1, etc.
    std::string name_;  

  public:
    // get the values of the composite parameter parameters
    std::vector<double> getParameterValues() const;
    const GeneCountDependencyFunction* getFunction() const{return func_;}

    // Returns true if the parameter is ignored
    static bool isIgnored(NcompositeParameter* param){return param == 0;}

    // get the function of the composite parameter
    const GeneCountDependencyFunction::FunctionType getFuncType() const {return func_->getName();}

    const std::string getName() const {return name_;}

    // get the number of parameters
    const size_t getSize() const {return params_.size();}

    // Returns the value of the overall independent/dependent rate on the number of chromosomes (for any function)
    double getRate(size_t state) const {return func_->getRate(params_, state);};

    // Converts a FunctionType enum variable to GeneCountDependencyFunction* object
    static GeneCountDependencyFunction* getDependencyFunction(GeneCountDependencyFunction::FunctionType funcType);


    static std::vector<std::string> getRelatedParameterNames(ParameterList &params, std::string pattern);
    NcompositeParameter(GeneCountDependencyFunction::FunctionType func, std::string paramName, vector<Parameter*> &params):
      params_(), func_(0), name_(paramName)
    {
      for (size_t i = 0; i < params.size(); i++){
        params_.push_back(params[i]);
      }
      //getParamUpdateFunction(func);
      func_ = getDependencyFunction(func);
    }


    virtual ~NcompositeParameter(){ delete func_;}



  protected:
    void setName(std::string name){name_ = name;}
    std::vector<Parameter*>& getParams(){return params_;}
    void setParams(std::vector<Parameter*> params){params_ = params;}
    friend GeneCountSubstitutionModel;  
};
  
  class GeneCountSubstitutionModel : public AbstractSubstitutionModel
  {
  public:
    enum paramType {LOSS, GAIN, INNOVATION, ELIMINATION};
    enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};

  private:
    NcompositeParameter* gain_;
    NcompositeParameter* loss_;
    NcompositeParameter* innovation_;
    NcompositeParameter* elimination_;
    int minState_;
    int maxState_;
    int countRange_;
    double firstNormQ_;
    mutable bool pijtCalledFromDeriv_;
    rootFreqType freqType_;
    static const std::map<int, std::string> newEventTypeToString;


  public:
    // Constructor
    GeneCountSubstitutionModel(const IntegerAlphabet* alpha, std::map<int, vector<double>> mapOfParamValues, rootFreqType freqType, ModelParameters* m);
    virtual ~GeneCountSubstitutionModel() {
      if (gain_ != 0){
        delete gain_;
      }
      if (loss_ != 0){
        delete loss_;
      }
      if (innovation_ != 0){
        delete innovation_;
      }
      if (elimination_ != 0){
        delete elimination_;
      }
    }
    GeneCountSubstitutionModel(const GeneCountSubstitutionModel& model):
      AbstractParameterAliasable(model),
      AbstractSubstitutionModel(model),
      gain_(0),
      loss_(0),
      innovation_(0),
      elimination_(0),
      minState_(model.minState_),
      maxState_(model.maxState_),
      countRange_(model.countRange_),
      freqType_(model.freqType_),
      firstNormQ_(model.firstNormQ_),
      pijtCalledFromDeriv_(model.pijtCalledFromDeriv_),
      vPowExp_(model.vPowExp_)
    {
      std::vector<NcompositeParameter**> newModelParams = {&gain_, &loss_, &innovation_, &elimination_};
      std::vector<NcompositeParameter*> originalModelParams = {model.gain_, model.loss_, model.innovation_, model.elimination_};
      for (size_t i = 0; i < originalModelParams.size(); i++){
        if (originalModelParams[i] == 0){
          continue;
        }
        vector<Parameter*> params = originalModelParams[i]->getParams();
        vector<Parameter*> newParams;
        for (size_t j = 0; j < params.size(); j++){
          std::string name = params[j]->getName();
          string noPrefixName = std::regex_replace(name, std::regex("GeneCount."), "");
          newParams.push_back(&(getParameter_(noPrefixName)));
        }

        *(newModelParams[i]) = new NcompositeParameter(originalModelParams[i]->getFuncType(), originalModelParams[i]->getName(), newParams);
        (*(newModelParams[i]))->func_->setDomainsIfNeeded(minState_, maxState_);

      }

    }

    GeneCountSubstitutionModel* clone() const override { return new GeneCountSubstitutionModel(*this); }
    std::string getName() const override { return "GeneCountSubstitutionModel"; }
    static int getParamIndexByName(string name);
    const Matrix<double>& getPij_t    (double d) const;
    const Matrix<double>& getdPij_dt  (double d) const;
    const Matrix<double>& getd2Pij_dt2(double d) const;
  
  protected:
    mutable std::vector< RowMatrix<double> > vPowExp_;
    std::vector<Parameter*> createRateParameter(GeneCountDependencyFunction::FunctionType &func, std::string paramName, vector<double> &vectorOfValues);

    // Creates, sets and adds all relevant substitution model parameters
    void updateParameters(std::map<int, vector<double>> mapOfParamValues, std::vector<int> rateChangeType);
    // Fills the generator matrix with parameter values
    void updateMatrices();

    void updateQWithGain(size_t i, size_t minChrNum);
    void updateQWithLoss(size_t i, size_t minChrNum);
    void updateQWithInnovation(size_t i);
    void updateQWithElimination(size_t i);
    double getFirstNorm() const;
    void calculateExp_Qt(size_t pow, double* s, double v) const;
    bool checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const;
    void updateEigenMatrices();
  };
}

#endif // GENECOUNT_SUBSTITUTION_MODEL_H