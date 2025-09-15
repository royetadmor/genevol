#include "GeneCountDependencyFunction.h"


using namespace std;
using namespace bpp;


/************************/
/*******CONSTANT*********/
/************************/
double ConstantDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue();
}

/************************/
/********LINEAR**********/
/************************/
double LinearDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  double func_res = params[0]->getValue() + ((double)(state-1)*params[1]->getValue());
  if (func_res < 0){
    return 0;
  }
  return func_res;

}

void LinearDependencyFunction::updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, int maxChrNum){
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

void LinearDependencyFunction::updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound){
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramName)))->getConstraint());
  interval->setLowerBound(lowerBound, interval->strictLowerBound());

}

void LinearDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
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

void LinearDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
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

/************************/
/*******LINEARBDD********/
/************************/
double LinearBDDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue() * (double)state;
}

/************************/
/*******EXPONENTIAL******/
/************************/
double ExponentailDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return params[0]->getValue() * std::exp((double)(state-1)*params[1]->getValue());
}

void ExponentailDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);

}

void ExponentailDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index > 1){
    throw Exception("ExponentailDependencyFunction::getAbsoluteBounds(): Too many parameters!!!");
  }
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam/(maxChrNumber-1);

  }
 
}

/************************/
/*******POLYNOMIAL*******/
/************************/
double PolynomialDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  return (params[0]->getValue()) * pow((double)(state) + params[1]->getValue(), params[2]->getValue());

}

void PolynomialDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);
}

void PolynomialDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  if (index == 0){
      *lowerBound = 0;
      *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = -domainMin_;
    *upperBound = upperBoundOfRateParam;
  }else if(index == 2){
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam;

  }else{
    throw Exception("PolynomialDependencyFunction::getAbsoluteBounds: index out of bounds!!");
    
  }

}

/************************/
/*******LOGNORMAL********/
/************************/
double LognormalDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  auto rangeFactor = params[0]->getValue();
  double scalingFactor = (double)domainMax_/logNormalDomainFactor;
  double transformedState = (double)state/scalingFactor;
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  auto eq_part_1 = 1/((double)(transformedState)*sigma*sqrt(2 * pi));
  auto eq_part_2 = std::exp(-(pow(log(transformedState)-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor *eq_part_1 * eq_part_2;

}

void LognormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);
}

void LognormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  *lowerBound = lowerBoundOfRateParam;
  if (index == 0){  // for the range parameter   
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // mu
    *upperBound = upperBoundLinearRateParam;
  }else if (index == 2){  // sigma
    *upperBound = upperBoundLinearRateParam*2;

  }else{
    throw Exception("LognormalDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
    
  }
}

/************************/
/******REVSIGMOID********/
/************************/
double RevSigmoidDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  // p1 is the range parameter
  auto p1 = params[0]->getValue();
  // p2 is the exponent multiplier parameter
  auto p2 = params[1]->getValue();
  // p3 is the shift parameter (should manipulate the cut of the reverse sigmoid tail)
  auto p3 = params[2]->getValue();
  // f(x) = p1* (e^-p2(x-p3)/(1+(e^-p2(x-p3))))
  auto x = static_cast<double>(state);
  return p1/(1+(std::exp(p2-(p3*x))));// (std::exp(-p2*(x-p3))/(1+(std::exp(-p2*(x-p3)))));

}

void RevSigmoidDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  
  if (index == 0){  // for the range parameter   
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // for the exponent parameter
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = (double)(domainMax_-domainMin_+1);
  }else if (index == 2){  // the shift parameter
    *lowerBound = -revSigmoidExpRateParam;
    *upperBound = lowerBoundOfRateParam;

  }else{
    throw Exception("RevSigmoidDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
    
  }
}

void RevSigmoidDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);
}

/************************/
/******LOGITNORMAL*******/
/************************/
double LogitnormalDependencyFunction::getRate(std::vector<Parameter*> params, size_t state) const{
  auto rangeFactor = params[0]->getValue();
  double transformedState = ((double)state-(double)domainMin_+1)/((double)domainMax_-(double)domainMin_+2);
  auto mu = params[1]->getValue();
  auto sigma = params[2]->getValue();
  double pi = 2 * acos(0.0);
  double expr_1 = 1/(sigma*sqrt(2 * pi));
  double expr_2 = 1/(transformedState*(1-transformedState));
  double logit_expr = log(transformedState/(1-transformedState));
  double expr_3 = std::exp(-(pow(logit_expr-mu, 2)/(2*pow(sigma, 2))));
  return rangeFactor * expr_1 * expr_2 * expr_3;

}

void LogitnormalDependencyFunction::getBoundsForInitialParams(size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber){
  getAbsoluteBounds(index, lowerBound, upperBound, maxChrNumber);
}

void LogitnormalDependencyFunction::getAbsoluteBounds(size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  *lowerBound = lowerBoundOfRateParam;
  if (index == 0){  // for the range parameter   
    *upperBound = upperBoundOfRateParam;
  }else if (index == 1){ // mu
    *upperBound = upperBoundLinearRateParam;
  }else if (index == 2){  // sigma
    *upperBound = upperBoundLinearRateParam*2;
  }else{
    throw Exception("LogitnormalDependencyFunction::getAbsoluteBounds(): index out of bounds!!");
  }
}

