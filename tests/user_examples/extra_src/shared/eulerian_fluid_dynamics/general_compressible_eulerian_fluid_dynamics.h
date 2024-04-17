/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2023 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	general_compressible_fluid_dynamics.h
 * @brief 	Here, we define the general compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */

#ifndef GENERAL_COMPRESSIBLE_EULERIAN_FLUID_DYNAMICS_H
#define GENERAL_COMPRESSIBLE_EULERIAN_FLUID_DYNAMICS_H

#include "compressible_fluid.h"
#include "fluid_body.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
  /**
  * @struct CompressibleFluidState
  * @brief  Struct for stored states of Riemann solver in compressible flow.
  */
  struct CompressibleFluidState : FluidState
  {
      Real& E_;
      CompressibleFluidState(Real& rho, Vecd& vel, Real& p, Real& E)
          : FluidState(rho, vel, p), E_(E) {};

  };

  /**
  * @struct DissipationState
  * @brief struct for stored dissipation state of Riemann solver for eulerian fluid dynamics.
  */
  struct DissipationStateCE
  {
      Matd momentum_dissipation_;
      Vecd density_dissipation_;
      Vecd energy_dissipation_;


      DissipationStateCE(Matd momentum_dissipation = Matd::Zero(), Vecd density_dissipation = Vecd::Zero(),
          Vecd energy_dissipation = Vecd::Zero()) : momentum_dissipation_(momentum_dissipation), 
           density_dissipation_(density_dissipation), energy_dissipation_(energy_dissipation) {};
  };

  /**
  * @struct NoRiemannSolverCE
  * @brief
  */
  class NoRiemannSolverCE
  {
  public:
      NoRiemannSolverCE(CompressibleFluid& compressible_fluid_i, CompressibleFluid& compressible_fluid_j) : 
          compressible_fluid_i_(compressible_fluid_i), compressible_fluid_j_(compressible_fluid_j) {};
      Vec2d getBoundingWaveSpeeds(const CompressibleFluidState& state_i, const CompressibleFluidState& state_j, const Vecd e_ij);
      DissipationStateCE getDissipationState(const CompressibleFluidState& state_i, const CompressibleFluidState& state_j, const Vecd e_ij);

  protected:
      CompressibleFluid & compressible_fluid_i_, & compressible_fluid_j_;
  };

  /**
  * @ struct HLLERiemannSolverCE
  * @ brief
  */
  class HLLERiemannSolverCE : public NoRiemannSolverCE
  {
  public:
      HLLERiemannSolverCE(CompressibleFluid& compressible_fluid_i, CompressibleFluid& compressible_fluid_j) : 
          NoRiemannSolverCE(compressible_fluid_i, compressible_fluid_j) {};
      Vec2d getBoundingWaveSpeeds(const CompressibleFluidState& state_i, const CompressibleFluidState& state_j, const Vecd e_ij);
      DissipationStateCE getDissipationState(const CompressibleFluidState& state_i, const CompressibleFluidState& state_j, const Vecd e_ij);
  };

  /*
   * @class EulerianCEViscousAccelerationInnerCE
   * @brief the viscosity force induced acceleration in Eulerian method
   */
  class EulerianCEViscousAccelerationInner : public fluid_dynamics::ViscousAccelerationInner
  {
  public:
      explicit EulerianCEViscousAccelerationInner(BaseInnerRelation& inner_relation);
      virtual ~EulerianCEViscousAccelerationInner() {};
      void interaction(size_t index_i, Real dt = 0.0);
      StdLargeVec<Real>& dE_dt_prior_;
      StdLargeVec<Vecd>& dmom_dt_prior_;

  };

  //----------------------------------------------------------------------
  //	Relaxation definition
  //----------------------------------------------------------------------
  /**
   * @class EulerianCompressibleTimeStepInitialization
   * @brief initialize a time step for a body.
   * including initialize particle acceleration
   * induced by viscous, gravity and other forces,
   * set the number of ghost particles into zero.
   */
  class EulerianCETimeStepInitialization : public TimeStepInitialization
  {
  public:
      EulerianCETimeStepInitialization(SPHBody& sph_body, SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd::Zero()));
      virtual ~EulerianCETimeStepInitialization() {};
      void update(size_t index_i, Real dt = 0.0);

  protected:
      StdLargeVec<Real>& rho_;
      StdLargeVec<Vecd>& pos_, & vel_;
      StdLargeVec<Vecd>& dmom_dt_prior_;
      StdLargeVec<Real>& dE_dt_prior_;
  };

  /**
 * @class EulerianCEAcousticTimeStepSize
 * @brief Computing the acoustic time step size
 */
  class EulerianCEAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
  {
  protected:
      StdLargeVec<Real>& rho_, & p_;
      StdLargeVec<Vecd>& vel_;
      Real smoothing_length_;

  public:
      explicit EulerianCEAcousticTimeStepSize(SPHBody& sph_body, Real acousticCFL = 0.6);
      virtual ~EulerianCEAcousticTimeStepSize() {};

      Real reduce(size_t index_i, Real dt = 0.0);
      virtual Real outputResult(Real reduced_value) override;
      CompressibleFluid compressible_fluid_;
      Real acousticCFL_;
  };

  /**
   * @class BaseIntegrationInGeneralCE
   * @brief Pure abstract base class for all fluid relaxation schemes in compressible flows
   */
  class BaseIntegrationInGeneralCE : public fluid_dynamics::BaseIntegration
  {
  public:
      explicit BaseIntegrationInGeneralCE(BaseInnerRelation& inner_relation);
      virtual ~BaseIntegrationInGeneralCE() {};

  protected:
      CompressibleFluid compressible_fluid_;
      StdLargeVec<Real>& Vol_, & E_, & dE_dt_, & dE_dt_prior_;
      StdLargeVec<Vecd>& mom_, & dmom_dt_, & dmom_dt_prior_;
  };

  /**
 * @class CEIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
  template <class RiemannSolverType>
  class CEIntegration1stHalf : public BaseIntegrationInGeneralCE
  {
  public:
      explicit CEIntegration1stHalf(BaseInnerRelation& inner_relation);
      virtual ~CEIntegration1stHalf() {};
      RiemannSolverType riemann_solver_;
      void interaction(size_t index_i, Real dt = 0.0);
      void update(size_t index_i, Real dt = 0.0);
  };
  using CEIntegration1stHalfNoRiemann = CEIntegration1stHalf<NoRiemannSolverCE>;
  using CEIntegration1stHalfHLLERiemann = CEIntegration1stHalf<HLLERiemannSolverCE>;

  /**
   * @class CEIntegration2ndHalf
   * @brief  Template density relaxation scheme in HLLC Riemann solver with and without limiter
   */
  template <class RiemannSolverType>
  class CEIntegration2ndHalf : public BaseIntegrationInGeneralCE
  {
  public:
      explicit CEIntegration2ndHalf(BaseInnerRelation& inner_relation);
      virtual ~CEIntegration2ndHalf() {};
      RiemannSolverType riemann_solver_;
      void interaction(size_t index_i, Real dt = 0.0);
      void update(size_t index_i, Real dt = 0.0);
  };
  using CEIntegration2ndHalfNoRiemann = CEIntegration2ndHalf<NoRiemannSolverCE>;
  using CEIntegration2ndHalfHLLERiemann = CEIntegration2ndHalf<HLLERiemannSolverCE>;

} // namespace SPH
#endif // GENERAL_COMPRESSIBLE_EULERIAN_FLUID_DYNAMICS_H