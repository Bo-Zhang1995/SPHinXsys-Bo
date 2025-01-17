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
 * @class BaseIntegration1stHalfCorrect
 */
template <class RiemannSolverType>
class BaseIntegration1stHalfCorrect : public BaseIntegration1stHalf<RiemannSolverType>
{
  public:
    explicit BaseIntegration1stHalfCorrect(BaseInnerRelation &inner_relation);
    virtual ~BaseIntegration1stHalfCorrect(){};

    using BaseIntegration1stHalf<RiemannSolverType>::BaseIntegration1stHalf;
    void initialization(size_t index_i, Real dt);
    void interaction(size_t index_i, Real dt);

  protected:
    StdLargeVec<Matd> p_B_;
    StdLargeVec<Matd> &B_;
};
using Integration1stHalfCorrect = BaseIntegration1stHalfCorrect<NoRiemannSolver>;
/** define the mostly used pressure relaxation scheme using Riemann solver */
using Integration1stHalfRiemannCorrect = BaseIntegration1stHalfCorrect<AcousticRiemannSolver>;

/**
 * @class BaseIntegration1stHalfConsistCorrect
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

} // namespace fluid_dynamics
} // namespace SPH
#endif // FLUID_DYNAMICS_INNER_CORRECTION_H