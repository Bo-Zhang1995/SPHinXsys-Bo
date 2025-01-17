/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	eulerian_fluid_dynamics.h
 * @brief 	Here, we define the common weakly compressible eulerian classes for fluid dynamics.
 * @author	Zhentong Wang and Xiangyu Hu
 */
#ifndef EULERIAN_FLUID_DYNAMICS_H
#define EULERIAN_FLUID_DYNAMICS_H

#include "fluid_body.h"
#include "fluid_dynamics_complex.h"
#include "fluid_dynamics_inner.h"
#include "general_dynamics.h"
#include "riemann_solver.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @struct EulerianNoRiemannSolver
 * @brief  Central difference scheme without Riemann flux.
 */
class EulerianNoRiemannSolver
{
  public:
    EulerianNoRiemannSolver(Fluid &fluid_i, Fluid &fluid_j)
        : fluid_i_(fluid_i), fluid_j_(fluid_j), 
        rho0_i_(fluid_i.ReferenceDensity()), rho0_j_(fluid_j.ReferenceDensity()),
        c0_i_(fluid_i.ReferenceSoundSpeed()), c0_j_(fluid_j.ReferenceSoundSpeed()),
        rho0c0_i_(rho0_i_ * c0_i_), rho0c0_j_(rho0_j_ * c0_j_),
        inv_rho0c0_sum_(1.0 / (rho0c0_i_ + rho0c0_j_)){};
    virtual FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij);
    Real DissipativePJump(const Real& u_jump);
    Real DissipativeUJump(const Real& p_jump, const Real& u_jump);

  protected:
    Fluid &fluid_i_, &fluid_j_;
    Real rho0_i_, rho0_j_;
    Real c0_i_, c0_j_;
    Real rho0c0_i_, rho0c0_j_, inv_rho0c0_sum_;
};

/**
 * @struct EulerianAcousticRiemannSolver
 * @brief  Acoustic RiemannSolver for Eulerian weakly-compressible flow.
 */
class EulerianAcousticRiemannSolver : public EulerianNoRiemannSolver
{
  public:
    EulerianAcousticRiemannSolver(Fluid &fluid_i, Fluid &fluid_j)
        : EulerianNoRiemannSolver(fluid_i, fluid_j),
        inv_rho0c0_ave_(2.0 * inv_rho0c0_sum_), 
        rho0c0_geo_ave_(2.0 * rho0c0_i_ * rho0c0_j_ * inv_rho0c0_sum_),
        inv_c_ave_(0.5 * (rho0_i_ + rho0_j_) * inv_rho0c0_ave_){};
    FluidStarState getInterfaceState(const FluidState &state_i, const FluidState &state_j, const Vecd &e_ij) override;
    Real DissipativePJump(const Real& u_jump);
    Real DissipativeUJump(const Real& p_jump, const Real& u_jump);

  protected:
    Real inv_rho0c0_ave_, rho0c0_geo_ave_;
    Real inv_c_ave_;
};

/**
 * @class EulerianIntegration1stHalf
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as template variable
 */
template <class RiemannSolverType>
class EulerianIntegration1stHalf : public BaseIntegration
{
  public:
    explicit EulerianIntegration1stHalf(BaseInnerRelation &inner_relation);
    virtual ~EulerianIntegration1stHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Vecd> &acc_prior_;
    StdLargeVec<Vecd> mom_, dmom_dt_;
};
/** define the mostly used pressure relaxation scheme using Riemann solver */
using EulerianIntegration1stHalfAcousticRiemann = EulerianIntegration1stHalf<EulerianAcousticRiemannSolver>;
using EulerianIntegration1stHalfNoRiemann = EulerianIntegration1stHalf<EulerianNoRiemannSolver>;

/**
 * @class EulerianIntegration1stHalfWithWall
 * @brief  template class pressure relaxation scheme with wall boundary
 */
template <class EulerianIntegration1stHalfType>
class EulerianIntegration1stHalfWithWall : public InteractionWithWall<EulerianIntegration1stHalfType>
{
  public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    EulerianIntegration1stHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration1stHalfType>(wall_contact_relation, base_body_relation){};
    explicit EulerianIntegration1stHalfWithWall(ComplexRelation &fluid_wall_relation)
        : EulerianIntegration1stHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~EulerianIntegration1stHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using EulerianIntegration1stHalfAcousticRiemannWithWall = EulerianIntegration1stHalfWithWall<EulerianIntegration1stHalfAcousticRiemann>;
using EulerianIntegration1stHalfNoRiemannWithWall = EulerianIntegration1stHalfWithWall<EulerianIntegration1stHalfNoRiemann>;

/**
 * @class EulerianIntegration2ndHalf
 * @brief  Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalf : public BaseIntegration
{
  public:
    explicit EulerianIntegration2ndHalf(BaseInnerRelation &inner_relation);
    virtual ~EulerianIntegration2ndHalf(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfAcousticRiemann = EulerianIntegration2ndHalf<EulerianAcousticRiemannSolver>;
using EulerianIntegration2ndHalfNoRiemann = EulerianIntegration2ndHalf<EulerianNoRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalfWithWall
 * @brief template density relaxation scheme with using  Riemann solver.
 */
template <class EulerianIntegration2ndHalfType>
class EulerianIntegration2ndHalfWithWall : public InteractionWithWall<EulerianIntegration2ndHalfType>
{
  public:
    // template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    EulerianIntegration2ndHalfWithWall(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration2ndHalfType>(wall_contact_relation, base_body_relation){};
    explicit EulerianIntegration2ndHalfWithWall(ComplexRelation &fluid_wall_relation)
        : EulerianIntegration2ndHalfWithWall(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~EulerianIntegration2ndHalfWithWall(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using EulerianIntegration2ndHalfAcousticRiemannWithWall = EulerianIntegration2ndHalfWithWall<EulerianIntegration2ndHalfAcousticRiemann>;
using EulerianIntegration2ndHalfNoRiemannWithWall = EulerianIntegration2ndHalfWithWall<EulerianIntegration2ndHalfNoRiemann>;

/*
 * @class EulerianIntegration1stHalfConsistency  
 * @brief Template class for pressure relaxation scheme with the Riemann solver
 * as temperate variable
 */
template <class RiemannSolverType>
class EulerianIntegration1stHalfConsistency : public BaseIntegration
{
  public:
    explicit EulerianIntegration1stHalfConsistency(BaseInnerRelation &inner_relation);
    virtual ~EulerianIntegration1stHalfConsistency(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
    StdLargeVec<Vecd> &acc_prior_;
    StdLargeVec<Vecd> mom_, dmom_dt_;
};
using EulerianIntegration1stHalfAcousticRiemannConsistency = EulerianIntegration1stHalfConsistency<EulerianAcousticRiemannSolver>;
using EulerianIntegration1stHalfNoRiemannConsistency = EulerianIntegration1stHalfConsistency<EulerianNoRiemannSolver>;

/**
 * @class EulerianIntegration1stHalfWithWallConsistency
 * @breif template class pressure relaxation scheme with wall boundary
 */
template <class EulerianIntegration1stHalfType>
class EulerianIntegration1stHalfWithWallConsistency : public InteractionWithWall<EulerianIntegration1stHalfType>
{
  public:
    template <class BaseBodyRelationType>
    EulerianIntegration1stHalfWithWallConsistency(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration1stHalfType>(wall_contact_relation, base_body_relation){};
    explicit EulerianIntegration1stHalfWithWallConsistency(ComplexRelation &fluid_wall_relation)
        : EulerianIntegration1stHalfWithWallConsistency(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~EulerianIntegration1stHalfWithWallConsistency(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using EulerianIntegration1stHalfAcousticRiemannWithWallConsistency =
    EulerianIntegration1stHalfWithWallConsistency<EulerianIntegration1stHalfAcousticRiemannConsistency>;
using EulerianIntegration1stHalfNoRiemannWithWallConsistency =
    EulerianIntegration1stHalfWithWallConsistency<EulerianIntegration1stHalfNoRiemannConsistency>;

/**
 * @class EulerianIntegration2ndHalfConsistency
 * @breif Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class EulerianIntegration2ndHalfConsistency : public BaseIntegration
{
  public:
    explicit EulerianIntegration2ndHalfConsistency(BaseInnerRelation &inner_relation);
    virtual ~EulerianIntegration2ndHalfConsistency(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using EulerianIntegration2ndHalfAcousticRiemannConsistency = EulerianIntegration2ndHalfConsistency<EulerianAcousticRiemannSolver>;
using EulerianIntegration2ndHalfNoRiemannConsistency = EulerianIntegration2ndHalfConsistency<EulerianNoRiemannSolver>;

/**
 * @class EulerianIntegration2ndHalfWithWallConsistency
 * @brief template density relaxation scheme with different Riemann solver
 */
template <class EulerianIntegration2ndHalfType>
class EulerianIntegration2ndHalfWithWallConsistency : public InteractionWithWall<EulerianIntegration2ndHalfType>
{
  public:
    //template for different combination of constructing body relations
    template <class BaseBodyRelationType>
    EulerianIntegration2ndHalfWithWallConsistency(BaseContactRelation &wall_contact_relation, BaseBodyRelationType &base_body_relation)
        : InteractionWithWall<EulerianIntegration2ndHalfType>(wall_contact_relation, base_body_relation){};
    explicit EulerianIntegration2ndHalfWithWallConsistency(ComplexRelation &fluid_wall_relation)
        : EulerianIntegration2ndHalfWithWallConsistency(fluid_wall_relation.getContactRelation(), fluid_wall_relation.getInnerRelation()){};
    virtual ~EulerianIntegration2ndHalfWithWallConsistency(){};
    void interaction(size_t index_i, Real dt = 0.0);
};
using EulerianIntegration2ndHalfAcousticRiemannWithWallConsistency =
    EulerianIntegration2ndHalfWithWallConsistency<EulerianIntegration2ndHalfAcousticRiemannConsistency>;
using EulerianIntegration2ndHalfNoRiemannWithWallConsistency =
    EulerianIntegration2ndHalfWithWallConsistency<EulerianIntegration2ndHalfNoRiemannConsistency>;

/**
 * @class SmearedSurfaceIndication
 * @brief Indication of the particles which are within cut-off radius of surface particles.
 */
class SmearedSurfaceIndication : public LocalDynamics, public FluidDataInner
{
  public:
    explicit SmearedSurfaceIndication(BaseInnerRelation &inner_relation);
    virtual ~SmearedSurfaceIndication(){};

    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<int> &indicator_;
    StdLargeVec<int> &smeared_surface_;
};

/**
 * @class NonReflectiveBoundaryCorrection
 * @brief Implement Eulerian non-reflective boundary condition at free surface particles.
 */
class NonReflectiveBoundaryCorrection : public LocalDynamics, public DataDelegateInner<BaseParticles>
{
  public:
    NonReflectiveBoundaryCorrection(BaseInnerRelation &inner_relation);
    virtual ~NonReflectiveBoundaryCorrection(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Fluid &fluid_;
    Real rho_farfield_, sound_speed_;
    Vecd vel_farfield_;
    StdLargeVec<Real> &rho_, &p_, &Vol_;
    StdLargeVec<Vecd> &vel_, &mom_, &pos_;
    StdLargeVec<Real> inner_weight_summation_, rho_average_, vel_normal_average_;
    StdLargeVec<Vecd> vel_tangential_average_, vel_average_;
    StdLargeVec<int> &indicator_, smeared_surface_;
    StdLargeVec<Vecd> &n_;
};
} // namespace fluid_dynamics
} // namespace SPH
#endif // EULERIAN_FLUID_DYNAMICS_H