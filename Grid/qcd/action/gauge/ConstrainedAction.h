/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ConstrainedAction.h

Copyright (C) 2026

Author: Ryan Hill <Ryan.Hill@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */
#pragma once

#include <cassert>
#include <sstream>

#include <Grid/qcd/action/ActionBase.h>

/*! \file
 *
 */
/// \cond DO_NOT_DOCUMENT
NAMESPACE_BEGIN(Grid);
/// \endcond

/*! @brief The parameters controlling the constrained action
 */
struct ConstrainedActionParameters
{
  RealD a;     ///< @brief LLR parameter - tuned to \f$d(ln[\rho(S)]) / dS\f$ at S0
  RealD S0;    ///< @brief Constrained action centre
  RealD sigma; ///< @brief Gaussian width
};

template <class WrappedAction>
class ConstrainedAction : public Action<typename WrappedAction::GaugeField>
/*! @brief Gaussian-constrained action for the LLR algorithm.
 *
 * Implements the constrained action
 *   \f$\beta S[U] \rightarrow aS[U] + (S[U]-S_0)^2/(2\sigma^2)\f$
 * given a gauge action S.
 * Assumes that \f$\beta\f$ for the unconstrained action has been set to 1
 * in order to perform the direct replacement of \f$S[U]\f$.
 *
 * Expects a single template parameter for the related unconstrained action.
 * @param WrappedAction: The unconstrained action type.
 */
{
public:
  using GaugeField = typename WrappedAction::GaugeField;

  using Action<GaugeField>::S;
  using Action<GaugeField>::Sinitial;
  using Action<GaugeField>::deriv;
  using Action<GaugeField>::refresh;

  /*! @brief Construct a constrained action.
   * @param[in] wrapped: the unconstrained action
   * @param[in] parameters: the constrained action parameter structure
   */
  ConstrainedAction(WrappedAction &wrapped, ConstrainedActionParameters parameters)
      : wrapped_(wrapped), parameters_(parameters)
  {
    this->is_smeared = wrapped_.is_smeared;
  }

  /*! @brief Delegate this to the unconstrained action.
   *
   * Gauge fields do not have pseudofermions, so this is a no-op
   */
  virtual void refresh(const GaugeField &U, GridSerialRNG &sRNG, GridParallelRNG &pRNG)
  {
    wrapped_.refresh(U, sRNG, pRNG);
  }

  /*! @brief The constrained gauge action itself
   * @param[in] U: The gauge field on which to compute the action.
   * @returns The value of the constrained action
   */
  virtual RealD S(const GaugeField &U)
  {
    return constrained_value(wrapped_.S(U));
  }

  /*! @brief The constrained value of the gauge action at the start of the trajectory.
   * @param[in] U: The gauge field on which to compute the action.
   * @returns The constrained value of the initial action
   */
  virtual RealD Sinitial(const GaugeField &U)
  {
    return constrained_value(wrapped_.Sinitial(U));
  }

  /*! @brief The related unconstrained gauge action
   * @param[in] U: The gauge field on which to compute the action.
   * @returns The value of the unconstrained action
   */
  RealD Sunconstrained(const GaugeField &U)
  {
    return wrapped_.S(U);
  }

  /*! @brief The derivative of the constrained gauge action
   *
   * This is the derivative of the unconstrained action
   * scaled by the appropriate a-dependent factor.
   *
   * @param[in] U: The gauge field on which to compute the derivative
   * @param[out] force: Output field into which to write the derivative
   */
  virtual void deriv(const GaugeField &U, GaugeField &force)
  {
    RealD base_action = wrapped_.S(U);
    wrapped_.deriv(U, force);

    RealD scale = parameters_.a + (base_action - parameters_.S0) / (parameters_.sigma * parameters_.sigma);
    force *= scale;
  }
  
  /*! @brief The constrained action name 
   * @returns The name of the constrained action
   */
  virtual std::string action_name()
  {
    return "ConstrainedAction<" + wrapped_.action_name() + ">";
  }

  /*! @brief A logger for the constrained action parameters. */
  virtual std::string LogParameters()
  {
    std::stringstream sstream;
    sstream << GridLogMessage << "[" << action_name() << "] a:     " << parameters_.a << std::endl;
    sstream << GridLogMessage << "[" << action_name() << "] S0:    " << parameters_.S0 << std::endl;
    sstream << GridLogMessage << "[" << action_name() << "] sigma: " << parameters_.sigma << std::endl;
    sstream << wrapped_.LogParameters();
    return sstream.str();
  }

  /*! @brief The constrained action parameters 
   * @returns The constrained action parameters structure
   */
  const ConstrainedActionParameters &parameters() const
  {
    return parameters_;
  }

  /*! @brief Set the constrained action parameters
   * @param[in] parameters: the constrained action parameter structure
  */
  void set_parameters(ConstrainedActionParameters parameters)
  {
    parameters_ = parameters;
  }

  /*! @brief Set the constrained action parameter \f$a\f$
   * @param[in] a: the constrained action parameter \f$a\f$ value
   */
  void set_a(RealD a)
  {
    parameters_.a = a;
  }

  /*! @brief Set the constrained action parameter \f$S_0\f$
   * @param[in] S0: the constrained action parameter \f$S_0\f$ value
   */
  void set_S0(RealD S0)
  {
    parameters_.S0 = S0;
  }

  /*! @brief Set the constrained action parameter \f$\sigma\f$
   * @param[in] sigma: the constrained action parameter \f$\sigma\f$ value
   */
  void set_sigma(RealD sigma)
  {
    parameters_.sigma = sigma;
  }

private:
  WrappedAction &wrapped_; ///< @brief The unconstrained action object
  ConstrainedActionParameters parameters_; ///< @brief The constrained action parameters structure

  /*! @brief Calculate the constrained action value from the unconstrained one.
  * @param[in] base_action: the unconstrained action value
  * @returns The constrained action value
  */
  RealD constrained_value(RealD base_action) const
  {
    RealD displacement = base_action - parameters_.S0;
    return (parameters_.a * base_action) + displacement * displacement / (2.0 * parameters_.sigma * parameters_.sigma);
  }
};

/// \cond DO_NOT_DOCUMENT
NAMESPACE_END(Grid);
/// \endcond

