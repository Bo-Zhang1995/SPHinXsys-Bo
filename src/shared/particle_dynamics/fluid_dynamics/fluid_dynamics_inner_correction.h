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
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file  fluid_dynamics_inner_correction.h
 * @brief Here, we define the algorithm classes for fluid dynamics,
 *        in which correction matrix is used to increase the approximation
 *        of pressure gradient.
 * @author Yaru Ren and Xiangyu Hu
 */

#ifndef FLUID_DYNAMICS_INNER_CORRECTION_H
#define FLUID_DYNAMICS_INNER_CORRECTION_H

#include "fluid_dynamics_inner.h"

namespace SPH
{
namespace fluid_dynamics
{
/**
 * @class BaseIntegration1stHalfConsistency
 */
template <class RiemannSolverType>
class BaseIntegration1stHalfConsistency : public BaseIntegration1stHalf<RiemannSolverType>
{
public:
    explicit BaseIntegration1stHalfConsistency(BaseInnerRelation& inner_relation);
    virtual ~BaseIntegration1stHalfConsistency() {};

    using BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf;
    void initialization(size_t index_i, Real dt);
    void interaction(size_t index_i, Real dt);

protected:
    StdLargeVec<Matd>& B_;
};
using Integration1stHalfConsistency = BaseIntegration1stHalfConsistency<NoRiemannSolver>;
/** define the mostly used pressure relaxation scheme using Riemann solver */
using Integration1stHalfRiemannConsistency = BaseIntegration1stHalfConsistency<AcousticRiemannSolver>;

/**
 * @class BaseIntegration2ndHalfCorrect
 * @brief Template density relaxation scheme with different Riemann solver
 */
template <class RiemannSolverType>
class BaseIntegration2ndHalfCorrect : public BaseIntegration2ndHalf<RiemannSolverType>
{
  public:
    explicit BaseIntegration2ndHalfCorrect(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration2ndHalfCorrect() {};

    using BaseIntegration2ndHalf<RiemannSolverType>::BaseIntegration2ndHalf;
    void initialization(size_t index_i, Real dt = 0.0);
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Matd> &B_;
};
using Integration2ndHalfCorrect = BaseIntegration2ndHalfCorrect<NoRiemannSolver>;
using Integration2ndHalfRiemannCorrect = BaseIntegration2ndHalfCorrect<AcousticRiemannSolver>;
} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_INNER_CORRECTION_H