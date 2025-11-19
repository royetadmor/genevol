#include "GeneCountSubstitutionModel.h"



using namespace bpp;
using namespace std;

const std::map<int, std::string> GeneCountSubstitutionModel::eventTypeToString = {
  {GeneCountSubstitutionModel::paramType::LOSS, "loss"},
  {GeneCountSubstitutionModel::paramType::GAIN, "gain"},
  {GeneCountSubstitutionModel::paramType::INNOVATION, "innovation"},
  {GeneCountSubstitutionModel::paramType::ELIMINATION, "elimination"}
};

// Constructor
GeneCountSubstitutionModel::GeneCountSubstitutionModel(std::shared_ptr<bpp::GeneCountAlphabet>& alpha,
  std::map<int, vector<double>> mapOfParamValues,
  rootFreqType freqType,
  std::vector<int> rateChangeType,
  ModelParameters* m):
    AbstractParameterAliasable("GeneCount."),
    AbstractSubstitutionModel(alpha, std::make_shared<const CanonicalStateMap>(alpha, false), "GeneCount."),
    gain_(0),
    loss_(0),
    innovation_(0),
    elimination_(0),
    freqType_(freqType),
    pijtCalledFromDeriv_(false)
{
  // Set class params
  minState_ = m->minState_;
  maxState_ = m->maxState_;
  countRange_ = m->countRange_;
  
  // Initiate generator matrix
  updateParameters(mapOfParamValues, rateChangeType);
  computeFrequencies(false);//TODO: should we actually do this?
  isScalable_ = false;
  updateMatrices_();
}

GeneCountDependencyFunction* compositeParameter::getDependencyFunction(GeneCountDependencyFunction::FunctionType funcType){
  switch (funcType)
  {
  case GeneCountDependencyFunction::CONSTANT:
    return new ConstantDependencyFunction();
  case GeneCountDependencyFunction::LINEAR:
    return new LinearDependencyFunction();
  case GeneCountDependencyFunction::EXP:
    return new ExponentialDependencyFunction();
  case GeneCountDependencyFunction::LINEAR_BD:
    return new LinearBDDependencyFunction();
  case GeneCountDependencyFunction::LOGNORMAL:
    return new LognormalDependencyFunction();
  case GeneCountDependencyFunction::POLYNOMIAL:
    return new PolynomialDependencyFunction();
  case GeneCountDependencyFunction::REVERSE_SIGMOID:
    return new RevSigmoidDependencyFunction();
  case GeneCountDependencyFunction::LOGITNORMAL:
    return new LogitnormalDependencyFunction();
  default:
    throw Exception("compositeParameter::getDependencyFunction(): No such function!!");
  }
}

void GeneCountSubstitutionModel::updateParameters(std::map<int, vector<double>> mapOfParamValues, std::vector<int> rateChangeType) {
  for ( int i = 0; i < NUM_MODEL_PARAMETERS; i++ )
  {
    auto funcType = static_cast<GeneCountDependencyFunction::FunctionType>(rateChangeType[i]);
    auto paramName = eventTypeToString.at(i);
    auto rateVector = mapOfParamValues.at(i);
    auto rParam = createRateParameter(funcType, paramName, rateVector);
    if (!rParam.empty()) {
      switch (i)
      {
      case GeneCountSubstitutionModel::paramType::GAIN:
        gain_ = new compositeParameter(funcType, paramName, rParam);
        break;
      case GeneCountSubstitutionModel::paramType::LOSS:
        loss_ = new compositeParameter(funcType, paramName, rParam);
        break;
      case GeneCountSubstitutionModel::paramType::INNOVATION:
        innovation_ = new compositeParameter(funcType, paramName, rParam);
        break;
      case GeneCountSubstitutionModel::paramType::ELIMINATION:
        elimination_ = new compositeParameter(funcType, paramName, rParam);
        break;
      default:
        std::string err = "GeneCountSubstitutionModel: no parameter named " + paramName;
        throw Exception(err);
      }
    }
  }
}

std::vector<Parameter*> GeneCountSubstitutionModel::createRateParameter(GeneCountDependencyFunction::FunctionType &func, std::string paramName, vector<double> &vectorOfValues){
  std::vector<Parameter*> params;
  if (func == GeneCountDependencyFunction::FunctionType::IGNORE){
    return params;
  }
  for (size_t i = 0; i < vectorOfValues.size(); i++){
    auto paramValue = vectorOfValues[i];
    if (paramValue == IgnoreParam){
      throw Exception("GeneCountSubstitutionModel::createRateParameter(): Function is not supposed to be defined as a legal function name if the parameter should be ignored!");
    }
    double lowerBound;
    double upperBound;
    GeneCountDependencyFunction* functionOp = compositeParameter::getDependencyFunction(func);
    functionOp->setDomainsIfNeeded(minState_, maxState_);
    functionOp->getAbsoluteBounds(i, &lowerBound, &upperBound, maxState_);
    //TODO: might be relevant, keeping for now
    // if ((simulated_) && (lowerBound > paramValue)){
    //   lowerBound = paramValue - EPSILON;

    // }
    delete functionOp;

    //TODO: might be relevant, keeping for now
    // if (simulated_ && upperBound <= paramValue){
    //   upperBound = paramValue + EPSILON;
    // }
    std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBound, upperBound, false, true);
    Parameter* param = new Parameter("GeneCount."+ paramName + std::to_string(i), paramValue, interval);
    params.push_back(param);
    addParameter_(param);
  }
  return params;
}

void GeneCountSubstitutionModel::updateMatrices_() {
 //update model parameters
  matchParametersValues(getParameters());
  
  //update generator matrix
  size_t maxChrNum = (size_t)(maxState_);
  size_t minChrNum = (size_t)(minState_); 
  MatrixTools::fill(generator_, 0);

  // updating Q matrix
  for (size_t i = minChrNum; i < maxChrNum+1; i++){
      if (i == 0) {
        updateQWithInnovation(i);
        continue;
      }
      if (i == 1) {
        updateQWithElimination(i);
      }
      if (i+1 < maxChrNum+1){
        updateQWithGain(i, minChrNum);
      }
      if (i > 1 && i-1 >= minChrNum){
        updateQWithLoss(i, minChrNum);
      }            
  }
  setDiagonal();
  updateEigenMatrices();
  firstNormQ_ = getFirstNorm();
}

void GeneCountSubstitutionModel::updateQWithGain(size_t i, size_t minChrNum){
  if (gain_ == 0){
    return;
  }
  double gainRate = gain_->getRate(i);
  if (gainRate < 0){
    throw Exception ("GeneCountSubstitutionModel::updateQWithGain(): negative gain rate!");
  }
  generator_(i-minChrNum, i+1-minChrNum) += gainRate;

}
/*******************************************************************************/
void GeneCountSubstitutionModel::updateQWithLoss(size_t i, size_t minChrNum){
  if (loss_ == 0){
    return;
  }
  double lossRate = loss_->getRate(i);
  if (lossRate < 0){
    throw Exception ("GeneCountSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  
  generator_(i-minChrNum, i-1-minChrNum) += lossRate;


}
/*******************************************************************************/
void GeneCountSubstitutionModel::updateQWithInnovation(size_t i){
    if (innovation_ == 0){
    return;
  }
  double innovationRate = innovation_->getRate(i);
  if (innovationRate < 0){
    throw Exception ("GeneCountSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  generator_(0, 1) = innovationRate;

}

void GeneCountSubstitutionModel::updateQWithElimination(size_t i){
  if (elimination_ == 0){
    return;
  }
  double eliminationRate = elimination_->getRate(i);
  if (eliminationRate < 0){
    throw Exception ("GeneCountSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  generator_(1, 0) = eliminationRate;

}

double GeneCountSubstitutionModel::getFirstNorm() const{
  double norm = 0;
  for (size_t i = 0; i < size_; i++){
    for (size_t j = 0; j < size_; j++){
      norm += fabs(generator_(i,j));

    }
  }
  return norm;
}

const Matrix<double>& GeneCountSubstitutionModel::getPij_t(double t) const
{
  size_t minTaylorIterations = 5;
  if (t == 0)
  {
    MatrixTools::getId(size_, pijt_);
  }
  else
  {
    RowMatrix<double> pijt_temp;
    MatrixTools::getId(size_, pijt_temp);
    double s = 1.0;
    double v = rate_ * t;
    double norm = v * firstNormQ_;
    size_t m = 0;
    bool converged = false;
    while (norm > 0.5)
    {
      m += 1;
      v /= 2;
      norm /= 2;
    }
    for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
      calculateExp_Qt(iternum, &s, v);

      if (iternum > minTaylorIterations){
        converged = checkIfReachedConvergence(pijt_, pijt_temp);
        if (converged){
          break;
        }
      }
      MatrixTools::copy(pijt_, pijt_temp);
      if (iternum > 250){
        throw Exception("GeneCountSubstitutionModel: Taylor series did not reach convergence!");
        break;
      }
      if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
        RowMatrix<double> new_pow;
        MatrixTools::mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
        vPowExp_.push_back(new_pow);
      }

    }
    while (m > 0){  // recover the 2^m
      
      MatrixTools::mult(pijt_, pijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, pijt_);

      m--;
    }
  }
  if (!pijtCalledFromDeriv_){
    for (size_t i = 0; i < size_; i++){
      for (size_t j = 0; j < size_; j++){
        if (pijt_(i,j) < 0){
          pijt_(i,j) = NumConstants::VERY_TINY(); // trying to do it exactly as in ChromEvol. Maybe the "nan" problem will be solved
          //pijt_(i,j) = 0;
        }
        else if (pijt_(i, j) > 1){
          pijt_(i,j) = 1.0;
        }

      }
    }

  }
  return pijt_;
}

const Matrix<double>& GeneCountSubstitutionModel::getdPij_dt(double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  pijt = getPij_t(d);
  MatrixTools::mult(pijt, generator_, dpijt_);
  MatrixTools::scale(dpijt_, rate_);
  pijtCalledFromDeriv_ = false;
  return dpijt_;

}

/******************************************************************************/
const Matrix<double>& GeneCountSubstitutionModel::getd2Pij_dt2(double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  pijt = getPij_t(d);
  MatrixTools::mult(vPowExp_[2], pijt, d2pijt_);
  MatrixTools::scale(d2pijt_, rate_ * rate_);
  pijtCalledFromDeriv_ = false;
  return d2pijt_;
}

void GeneCountSubstitutionModel::calculateExp_Qt(size_t pow, double* s, double v) const{
  if (pow == 2){
    MatrixTools::getId(size_, pijt_);
    for (size_t i = 1; i <= pow; i++) {
      *s *= v / static_cast<double>(i);// the initial value of v is rt/(2^m)
      MatrixTools::add(pijt_, *s, vPowExp_[i]);
    }
  } else {
    *s *= v / static_cast<double>(pow);
    MatrixTools::add(pijt_, *s, vPowExp_[pow]);
  }
}

bool GeneCountSubstitutionModel::checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const{
    for (size_t i = 0; i < pijt.getNumberOfRows(); i++){
        for (size_t j = 0; j < pijt.getNumberOfColumns(); j++){
            double diff = fabs(pijt(i,j) - mt_prev(i,j));
            if (diff > TAYLOR_ITERATION_EPSILON){
                return false;
            }else if ((pijt(i,j) + TAYLOR_ITERATION_EPSILON < 0) || (pijt(i,j) > 1 + TAYLOR_ITERATION_EPSILON)){
              return false;
            }
        }
    }
    return true;
}

void GeneCountSubstitutionModel::updateEigenMatrices()
{
  // Compute eigen values and vectors:
  if (enableEigenDecomposition())
  {
    // Look for null lines (such as stop lines)
    // ie null diagonal elements

    size_t nbStop=0;
    size_t salph = getNumberOfStates();
    vector<bool> vnull(salph); // vector of the indices of lines with
                               // only zeros

    for (size_t i = 0; i < salph; i++)
    {
      if (abs(generator_(i, i)) < NumConstants::TINY())
      {
        nbStop++;
        vnull[i]=true;
      }
      else
        vnull[i]=false;
    }
        
    if (nbStop != 0)
    {
      size_t salphok=salph - nbStop;
      
      RowMatrix<double> gk(salphok, salphok);
      size_t gi = 0, gj = 0;

      for (size_t i = 0; i < salph; i++)
      {
        if (!vnull[i])
        {
          gj = 0;
          for (size_t j = 0; j < salph; j++)
          {
            if (!vnull[j])
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();

      for (size_t i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
        iEigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (size_t i = 0; i < salph; i++)
      {
        if (vnull[i])
        {
          gi++;
          for (size_t j = 0; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }

          rightEigenVectors_(i, salphok + gi - 1) = 1;
        }
        else
        {
          for (size_t j = 0; j < salphok; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi, j);
          }

          for (size_t j = salphok; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      rightEigenVectors_ = ev.getV();
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();
      nbStop = 0;
    }

    /// Now check inversion and diagonalization
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

      // is it diagonalizable ?
      isDiagonalizable_ = true;

      if (!dynamic_cast<ReversibleSubstitutionModelInterface*>(this))
      {
        for (auto& vi : iEigenValues_)
        {
          if (abs(vi) > NumConstants::TINY())
          {
            isDiagonalizable_ = false;
            break;
          }
        }
      }
      
      // looking for the vector of 0 eigenvalues

      std::vector<size_t> vNullEv;
      for (size_t i = 0; i< salph - nbStop; i++)
        if ((abs(eigenValues_[i]) < NumConstants::SMALL()) && (abs(iEigenValues_[i]) < NumConstants::SMALL()))
          vNullEv.push_back(i);
      

      // pb to find unique null eigenvalue      
      isNonSingular_=(vNullEv.size()==1);

      size_t nulleigen;
      
      double val;
      if (!isNonSingular_)
      {
        //look or check which non-stop right eigen vector elements are
        //equal.
        for (auto cnull : vNullEv)
        {
          size_t i = 0;
          while (vnull[i])
            i++;
          
          val = rightEigenVectors_(i, cnull);
          i++;
          
          while (i < salph)
          {
            if (!vnull[i])
            {
              if (abs(rightEigenVectors_(i, cnull) - val) > NumConstants::SMALL())
                break;
            }
            i++;
          }
          
          if (i >= salph)
          {
            isNonSingular_ = true;
            nulleigen=cnull;
            break;
          }
        }
      }
      else
        nulleigen=vNullEv[0];
      
      if (isNonSingular_)
      {
        eigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
        iEigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches


      }
      else
      {
        //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Unable to find eigenvector for eigenvalue 0. Taylor series used instead.");
        isDiagonalizable_ = false;
      }
    }
    // if rightEigenVectors_ is singular
    catch (ZeroDivisionException& e)
    {
      //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Singularity during diagonalization. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (vPowExp_.size() == 0){
      vPowExp_.resize(30);
    }
    MatrixTools::getId(salph, vPowExp_[0]);
    MatrixTools::Taylor(generator_, 30, vPowExp_);

  }

}

int GeneCountSubstitutionModel::getParamIndexByName(string name) {
  for (const auto& kv : eventTypeToString)
  {
      if (name.find(kv.second) != std::string::npos)
            return kv.first;
  }
  throw std::runtime_error("Value not found in map: " + name);
}
