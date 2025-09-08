#include "GeneCountDependencyFunction.h"


using namespace std;
using namespace bpp;


double NConstantDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue();
}

/**************************************************************************************/
double NLinearDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  double func_res = params[0]->getValue() + ((double)(state-1)*params[1]->getValue());
  if (func_res < 0){
    return 0;
  }
  return func_res;

}
/**************************************************************************************/
void NLinearDependencyFunction::updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum){
  if (index == 0){
    *lowerBound = std::max(lowerBoundOfRateParam, -params.getParameter(paramsNames[1]).getValue()*(maxChrNum-1));
    *upperBound = upperBoundOfRateParam;
    
  }else if (index == 1){
    *lowerBound = -params.getParameter(paramsNames[0]).getValue()/(maxChrNum-1);
    *upperBound = upperBoundLinearRateParam;   
  }else{
    throw Exception("LinearDependencyFunction::updateBounds(): index out of bounds!!");
    
  }
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
  interval->setLowerBound(*lowerBound, interval->strictLowerBound());
}
/**************************************************************************************/
void NLinearDependencyFunction::updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramName)))->getConstraint());
  interval->setLowerBound(lowerBound, interval->strictLowerBound());

}
/**************************************************************************************/
void NLinearDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;
    
  }else if (index == 1){
    *lowerBound = -paramValues[0]/(maxChrNumber-1);
    *upperBound = upperBoundLinearRateParam;       
  }else{
    throw Exception("LinearDependencyFunction::updateBounds(): index out of bounds!!");
    
  }
}
/**************************************************************************************/
void NLinearDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index == 0){
    *lowerBound = -upperBoundLinearRateParam*(maxChrNumber - 1);
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){
    *lowerBound = -upperBoundOfRateParam/(maxChrNumber -1);
    *upperBound = upperBoundLinearRateParam;
  }else{
    throw Exception("LinearDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
  }
}
/**************************************************************************************/