#ifndef EXTENDED_BRENT_OPTIMIZER_H
#define EXTENDED_BRENT_OPTIMIZER_H


#include <Bpp/Numeric/Function/AbstractOptimizer.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/NumTools.h>
#include <Bpp/Numeric/Function/OneDimensionOptimizationTools.h>

namespace bpp
{
class ExtendedBrentOptimizer :
  public AbstractOptimizer
{
public:
  class BODStopCondition :
    public AbstractOptimizationStopCondition
  {
public:
    BODStopCondition(ExtendedBrentOptimizer* bod) :
      AbstractOptimizationStopCondition(bod)
    {
      tolerance_ = bod->tol2;
      burnin_ = 3;
    }
    virtual ~BODStopCondition() {}

    BODStopCondition* clone() const { return new BODStopCondition(*this); }

public:
    bool isToleranceReached() const;
    double getCurrentTolerance() const;
  };

  enum Bracketing
  {
    BRACKET_OUTWARD = 0,
    BRACKET_INWARD = 1,
    BRACKET_SIMPLE = 2
  };

  friend class BODStopCondition;

protected:
  double a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
  double _xinf, _xsup;
  bool isInitialIntervalSet_;
  Bracketing bracketing_;

public:
  ExtendedBrentOptimizer(std::shared_ptr<FunctionInterface> function = nullptr);
  virtual ~ExtendedBrentOptimizer() {}

  ExtendedBrentOptimizer* clone() const { return new ExtendedBrentOptimizer(*this); }

public:
  /**
   * @name The Optimizer interface.
   *
   * @{
   */

  /**
   * @brief Initialize optimizer.
   *
   * Brent's algorithm needs 2 initial guesses, so you must call the
   * setInitialInterval() method first. This function actually performs:
   * <ul>
   * <li>Parameter list actualisation;</li>
   * <li>Initial bracketting;</li>
   * <li>Function evaluation count reseting.</li>
   * </ul>
   */
  double optimize(); // redefinition
  /** @} */

  void doInit(const ParameterList& params);

  double doStep();

  Bracket setSimpleBracketing(double a, double c, FunctionInterface& function, ParameterList parameters);

public:
  /**
   * @name Specific method
   *
   * @{
   */

  /**
   * @brief Set intial search interval.
   *
   * @param inf Minimum value.
   * @param sup Maximum value.
   */
  void setInitialInterval(double inf, double sup);
  /** @} */

  /**
   * @return 'true' if the initial interval has been correctly set.
   */
  bool isInitialIntervalSet() const { return isInitialIntervalSet_; }

  /**
   * @brief Get the brackeitng method
   */
  ExtendedBrentOptimizer::Bracketing getBracketing() const { return bracketing_; }

  /**
   * @brief Set the brackeitng method
   */
  void setBracketing(ExtendedBrentOptimizer::Bracketing bracketing)  { bracketing_ = bracketing; }

public:
  static double ZEPS;
};
} // end of namespace bpp.
#endif // EXTENDED_BRENT_OPTIMIZER_H