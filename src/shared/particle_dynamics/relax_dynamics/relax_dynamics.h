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
 * @file 	relax_dynamics.h
 * @brief 	This is the classes of particle relaxation in order to produce body fitted
 * 			initial particle distribution.
 * @author	Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef RELAX_DYNAMICS_H
#define RELAX_DYNAMICS_H

#include "all_body_relations.h"
#include "all_particle_dynamics.h"
#include "base_kernel.h"
#include "cell_linked_list.h"
#include "general_dynamics.h"

namespace SPH
{
class GeometryShape;
class LevelSetShape;

class PressureRelaxation
{
  public:
    PressureRelaxation(){};

    Real getBackgroundForce(Matd Bi, Matd Bj)
    {
        return 2.0;
    };
};

class CorrectionMatrixRelaxation
{
  public:
    CorrectionMatrixRelaxation(){};

    Matd getBackgroundForce(Matd Bi, Matd Bj)
    {
        return (Bi + Bj);
    };
};

namespace relax_dynamics
{
typedef DataDelegateSimple<BaseParticles> RelaxDataDelegateSimple;

typedef DataDelegateInner<BaseParticles> RelaxDataDelegateInner;

typedef DataDelegateComplex<BaseParticles, BaseParticles> RelaxDataDelegateComplex;

/**
 * @class GetTimeStepSizeSquare
 * @brief relaxation dynamics for particle initialization
 * computing the square of time step size
 */
class GetTimeStepSizeSquare : public LocalDynamicsReduce<Real, ReduceMax>,
                              public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Vecd> &acc_;
    Real h_ref_;

  public:
    explicit GetTimeStepSizeSquare(SPHBody &sph_body);
    virtual ~GetTimeStepSizeSquare(){};

    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value);
};

/**
 * @class RelaxationAccelerationInner
 * @brief simple algorithm for physics relaxation
 * without considering contact interaction.
 * this is usually used for solid like bodies
 */
template <typename RelaxationType>
class RelaxationAccelerationInner : public LocalDynamics, 
                                    public RelaxDataDelegateInner
{
  public:
    explicit RelaxationAccelerationInner(BaseInnerRelation &inner_relation);
    virtual ~RelaxationAccelerationInner(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            acceleration -= relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) * 
                            inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }
        acc_[index_i] = acceleration;
    };

  protected:
    RelaxationType relaxation_type;
    StdLargeVec<Vecd> &acc_, &pos_;
    StdLargeVec<Matd>& B_;
};

/**
 * @class RelaxationAccelerationInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
template <typename RelaxationType>
class RelaxationAccelerationInnerWithLevelSetCorrection : public RelaxationAccelerationInner<RelaxationType>
{
  public:
    explicit RelaxationAccelerationInnerWithLevelSetCorrection(BaseInnerRelation &inner_relation);
    virtual ~RelaxationAccelerationInnerWithLevelSetCorrection(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationInner<RelaxationType>::interaction(index_i, dt);
        Real overlap = level_set_shape_->computeKernelIntegral(this->pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));
       
        /* A scaling is adopted to handle the particle overlap. */
        this->acc_[index_i] -= this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) * 
                         level_set_shape_->computeKernelGradientIntegral(this->pos_[index_i],
                         sph_adaptation_->SmoothingLengthRatio(index_i));
    };

  protected:
    LevelSetShape *level_set_shape_;
    SPHAdaptation *sph_adaptation_;
};

/**
 * @class UpdateParticlePosition
 * @brief update the particle position for a time step
 */
class UpdateParticlePosition : public LocalDynamics,
                               public RelaxDataDelegateSimple
{
  protected:
    SPHAdaptation *sph_adaptation_;
    StdLargeVec<Vecd> &pos_, &acc_;

  public:
    explicit UpdateParticlePosition(SPHBody &sph_body);
    virtual ~UpdateParticlePosition(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class UpdateSmoothingLengthRatioByShape
 * @brief update the particle smoothing length ratio
 */
class UpdateSmoothingLengthRatioByShape : public LocalDynamics,
                                          public RelaxDataDelegateSimple
{
  protected:
    StdLargeVec<Real> &h_ratio_, &Vol_;
    StdLargeVec<Vecd> &pos_;
    Shape &target_shape_;
    ParticleRefinementByShape *particle_adaptation_;
    Real reference_spacing_;

  public:
    UpdateSmoothingLengthRatioByShape(SPHBody &sph_body, Shape &target_shape);
    explicit UpdateSmoothingLengthRatioByShape(SPHBody &sph_body);
    virtual ~UpdateSmoothingLengthRatioByShape(){};

    void update(size_t index_i, Real dt = 0.0);
};

/**
 * @class ShapeSurfaceBounding
 * @brief constrain surface particles by
 * map constrained particles to geometry face and
 * r = r + phi * norm (vector distance to face)
 */
class ShapeSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                             public RelaxDataDelegateSimple
{
  public:
    ShapeSurfaceBounding(NearShapeSurface &body_part);
    virtual ~ShapeSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    LevelSetShape *level_set_shape_;
    Real constrained_distance_;
};

/**
 * @class RelaxationStepInner
 * @brief carry out particle relaxation step of particles within the body
 */
template <typename RelaxationType = PressureRelaxation>
class RelaxationStepInner : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepInner(BaseInnerRelation &inner_relation,
                                 bool level_set_correction = false);
    virtual ~RelaxationStepInner(){};
    SimpleDynamics<ShapeSurfaceBounding> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    BaseInnerRelation &inner_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_inner_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
 * @class RelaxationAccelerationComplex
 * @brief compute relaxation acceleration while consider the present of contact bodies
 * with considering contact interaction
 * this is usually used for fluid like bodies //TODO: seems better called as Contact
 */
template <typename RelaxationType = PressureRelaxation>
class RelaxationAccelerationComplex : public LocalDynamics,
                                      public RelaxDataDelegateComplex
{
public:
    explicit RelaxationAccelerationComplex(ComplexRelation& complex_relation);
    virtual ~RelaxationAccelerationComplex() {};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        Vecd acceleration = Vecd::Zero();
        Neighborhood& inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            acceleration -= relaxation_type.getBackgroundForce(B_[index_i], B_[index_j]) *
                            inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
        }

        /** Contact interaction. */
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            StdLargeVec<Matd>& B_k = *(contact_B_[k]);
            Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                acceleration -= relaxation_type.getBackgroundForce(B_[index_i], B_k[index_j]) *
                                contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];
            }
        }

        acc_[index_i] = acceleration;
    };

protected:
    RelaxationType relaxation_type;
    StdLargeVec<Vecd>& acc_, & pos_;
    StdLargeVec<Matd>& B_;
    StdVec<StdLargeVec<Matd>*> contact_B_;
};

/**
 * @class RelaxationAccelerationComplexWithLevelSetCorrection
 * @brief compute relaxation acceleration while consider the present of contact bodies
 * with considering contact interaction
 * this is usually used for fluid like bodies
 * we constrain particles with a level-set correction function when the fluid boundary is not contacted with solid.
 */
template <typename RelaxationType = PressureRelaxation>
class RelaxationAccelerationComplexWithLevelSetCorrection : public RelaxationAccelerationComplex<RelaxationType>
{
  public:
    RelaxationAccelerationComplexWithLevelSetCorrection(
        ComplexRelation &complex_relation, const std::string &shape_name);
    virtual ~RelaxationAccelerationComplexWithLevelSetCorrection(){};

    inline void interaction(size_t index_i, Real dt = 0.0)
    {
        RelaxationAccelerationComplex<RelaxationType>::interaction(index_i, dt);
        Real overlap = level_set_shape_->computeKernelIntegral(this->pos_[index_i], 
                       sph_adaptation_->SmoothingLengthRatio(index_i));

        this->acc_[index_i] -= this->relaxation_type.getBackgroundForce(this->B_[index_i], this->B_[index_i]) *
                               level_set_shape_->computeKernelGradientIntegral(this->pos_[index_i],
                               sph_adaptation_->SmoothingLengthRatio(index_i));
    };

  protected:
    LevelSetShape *level_set_shape_;
    SPHAdaptation *sph_adaptation_;
};

/**
 * @class RelaxationStepComplex
 * @brief carry out particle relaxation step of particles within multi bodies
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationStepComplex : public BaseDynamics<void>
{
  public:
    explicit RelaxationStepComplex(ComplexRelation &complex_relation,
                                   const std::string &shape_name, bool level_set_correction = false);
    virtual ~RelaxationStepComplex(){};
    SimpleDynamics<ShapeSurfaceBounding> &SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

  protected:
    RealBody *real_body_;
    ComplexRelation &complex_relation_;
    NearShapeSurface near_shape_surface_;
    UniquePtr<BaseDynamics<void>> relaxation_acceleration_complex_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_square_;
    SimpleDynamics<UpdateParticlePosition> update_particle_position_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

template <typename ErrorDataType, typename ParameterADataType, typename ParameterCDataType>
struct ErrorAndParameters
{
    ErrorDataType error_;
    ParameterADataType a_;
    ParameterCDataType c_;

    ErrorAndParameters<ErrorDataType, ParameterADataType, ParameterCDataType>() :
        error_(ZeroData<ErrorDataType>::value),
        a_(ZeroData<ParameterADataType>::value),
        c_(ZeroData<ParameterCDataType>::value) {};
};

/**
 * @class RelaxationInnerImplicit
 * @brief carry out particle relaxation by position with implicit evolution.
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationInnerImplicit : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    explicit RelaxationInnerImplicit(BaseInnerRelation& inner_relation);
    virtual ~RelaxationInnerImplicit() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> implicit_residual_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    RelaxationType relaxation_type;
};

/**
 * @class RelaxationByCMImplicitInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationInnerWithLevelSetCorrectionImplicit : public RelaxationInnerImplicit<RelaxationType>
{
public:
    explicit RelaxationInnerWithLevelSetCorrectionImplicit(BaseInnerRelation& inner_relation);
    virtual ~RelaxationInnerWithLevelSetCorrectionImplicit() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepInnerImplicit
 * @brief carry out the particle relaxation evolution from first order consistency within the body
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationStepInnerImplicit : public BaseDynamics<void>
{
public:
    explicit RelaxationStepInnerImplicit(BaseInnerRelation& inner_relation,
                                         bool level_set_correction = false);
    virtual ~RelaxationStepInnerImplicit() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real time_step_size_;
    RealBody* real_body_;
    BaseInnerRelation& inner_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationInnerWithLevelSetCorrectionImplicit<RelaxationType>> relaxation_evolution_inner_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
 * @class RelaxationComplexImplicit
 * @brief carry out particle relaxation by position with implicit evolution.
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationComplexImplicit : public LocalDynamics, public RelaxDataDelegateComplex
{
public:
    explicit RelaxationComplexImplicit(ComplexRelation& complex_relation, const std::string& shape_name);
    virtual ~RelaxationComplexImplicit() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0);
    virtual void updateStates(size_t index_i, Real dt, const ErrorAndParameters<Vecd, Matd, Matd>& error_and_parameters);

    Kernel* kernel_;
    StdLargeVec<Real>& Vol_;
    StdVec<StdLargeVec<Real>*> contact_Vol_;
    StdLargeVec<Vecd>& pos_, & acc_;
    StdLargeVec<Matd>& B_;
    StdVec<StdLargeVec<Matd>*> contact_B_;
    StdLargeVec<Real> &implicit_residual_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    RelaxationType relaxation_type;
};

/**
 * @class RelaxationByCMImplicitInnerWithLevelSetCorrection
 * @brief we constrain particles to a level function representing the interface.
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationComplexWithLevelSetCorrectionImplicit : public RelaxationComplexImplicit<RelaxationType>
{
public:
    explicit RelaxationComplexWithLevelSetCorrectionImplicit(ComplexRelation& complex_relation, const std::string& shape_name);
    virtual ~RelaxationComplexWithLevelSetCorrectionImplicit() {};

protected:
    virtual ErrorAndParameters<Vecd, Matd, Matd> computeErrorAndParameters(size_t index_i, Real dt = 0.0) override;
};

/**
 * @class RelaxationStepComplexImplicit
 * @brief carry out the particle relaxation evolution from first order consistency within the body
 */
template <class RelaxationType = PressureRelaxation>
class RelaxationStepComplexImplicit : public BaseDynamics<void>
{
public:
    explicit RelaxationStepComplexImplicit(ComplexRelation& complex_relation, const std::string& shape_name,
                                           bool level_set_correction = false);
    virtual ~RelaxationStepComplexImplicit() {};
    SimpleDynamics<ShapeSurfaceBounding>& SurfaceBounding() { return surface_bounding_; };
    virtual void exec(Real dt = 0.0) override;

protected:
    Real time_step_size_;
    RealBody* real_body_;
    ComplexRelation& complex_relation_;
    NearShapeSurface near_shape_surface_;
    ReduceDynamics<GetTimeStepSizeSquare> get_time_step_;
    InteractionSplit<RelaxationComplexImplicit<RelaxationType>> relaxation_evolution_complex_;
    SimpleDynamics<ShapeSurfaceBounding> surface_bounding_;
};

/**
 * @class UpdateParticleKineticEnergy
 * @brief calculate the particle kinetic energy
 */
class UpdateParticleKineticEnergy : public LocalDynamics, public RelaxDataDelegateInner
{
public:
    UpdateParticleKineticEnergy(BaseInnerRelation& inner_relation);
    virtual ~UpdateParticleKineticEnergy() {};
    void update(size_t index_i, Real dt);

protected:
    StdLargeVec<Real>& mass_;
    StdLargeVec<Vecd>& acc_;
    StdLargeVec<Real> particle_kinetic_energy;
};

/**
 * @class CheckCorrectedZeroOrderConsistency
 * @brief calculate the corrected zero order consistency
 */
class CheckCorrectedZeroOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedZeroOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedZeroOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_zero_order_error_norm_;
    StdLargeVec<Vecd> corrected_zero_order_error_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class CheckCorrectedFirstOrderConsistency
 * @brief calculate the corrected first order consistency
 */
class CheckCorrectedFirstOrderConsistency : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckCorrectedFirstOrderConsistency(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckCorrectedFirstOrderConsistency() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> corrected_first_order_error_norm_;
    StdLargeVec<Matd> corrected_first_order_error_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class CheckL2NormError
 * @brief
 */
class CheckL2NormError : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckL2NormError(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckL2NormError() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& scalar_;
    StdLargeVec<Real> average_label_;
    StdLargeVec<Real> analytical_;
    StdLargeVec<Real> L2_NKGC_;
    StdLargeVec<Real> L2_SKGC_;
    StdLargeVec<Real> L2_CKGC_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Vecd>& pos_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
    Real constrained_distance_;
};

/**
 * @class CheckConsistencyRealization
 * @brief check the consistency of SPH conservative formulation.
 */
class CheckConsistencyRealization : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckConsistencyRealization() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& scalar_;
    StdLargeVec<Vecd>& vector_;
    StdLargeVec<Matd>& matrix_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Vecd> reproduce_scalar_gradient_;
    StdLargeVec<Real> reproduce_vector_gradient_;
    StdLargeVec<Vecd> reproduce_matrix_gradient_;
    StdLargeVec<Real> ACterm_norm_;
    StdLargeVec<Vecd> ACterm_;
    StdLargeVec<Real> ASterm_norm_;
    StdLargeVec<Vecd> ASterm_;
    StdLargeVec<Real> Cterm_norm_;
    StdLargeVec<Vecd> Cterm_;
    StdLargeVec<Real> NKGC_norm_;
    StdLargeVec<Vecd> NKGC_;
    StdLargeVec<Real> PR_ERROR_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class CheckReverseConsistencyRealization
 */
class CheckReverseConsistencyRealization : public LocalDynamics, public GeneralDataDelegateInner
{
public:
    CheckReverseConsistencyRealization(BaseInnerRelation& inner_relation, bool level_set_correction = false);
    virtual ~CheckReverseConsistencyRealization() {};
    void interaction(size_t index_i, Real dt);

protected:
    bool level_set_correction_;
    StdLargeVec<Real>& scalar_;
    StdLargeVec<Vecd>& pos_;
    StdLargeVec<Matd>& B_;
    StdLargeVec<Real> ABterm_norm_;
    StdLargeVec<Vecd> ABterm_;
    StdLargeVec<Real> ARterm_norm_;
    StdLargeVec<Vecd> ARterm_;
    StdLargeVec<Real> Bterm_norm_;
    StdLargeVec<Vecd> Bterm_;
    LevelSetShape* level_set_shape_;
    SPHAdaptation* sph_adaptation_;
};

/**
 * @class ShellMidSurfaceBounding
 * @brief constrain particles by constraining particles to mid-surface.
 * Note that level_set_refinement_ratio should be smaller than particle_spacing_ref_ / (0.05 * thickness_)
 * because if level_set_refinement_ratio > particle_spacing_ref_ / (0.05 * thickness_),
 * there will be no level set field.
 */
class ShellMidSurfaceBounding : public BaseLocalDynamics<BodyPartByCell>,
                                public RelaxDataDelegateInner
{
  public:
    ShellMidSurfaceBounding(NearShapeSurface &body_part, BaseInnerRelation &inner_relation);
    virtual ~ShellMidSurfaceBounding(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    StdLargeVec<Vecd> &pos_;
    Real constrained_distance_;
    Real particle_spacing_ref_;
    LevelSetShape *level_set_shape_;
};

/**
 * @class ShellNormalDirectionPrediction
 * @brief predict the normal direction of shell particles.
 */
class ShellNormalDirectionPrediction : public BaseDynamics<void>
{
    const Real convergence_criterion_;
    const Real consistency_criterion_;

    void predictNormalDirection();
    void correctNormalDirection();

  public:
    explicit ShellNormalDirectionPrediction(BaseInnerRelation &inner_relation,
                                            Real thickness, Real consistency_criterion = cos(Pi / 20.0));
    virtual ~ShellNormalDirectionPrediction(){};
    virtual void exec(Real dt = 0.0) override;

  protected:
    class NormalPrediction : public RelaxDataDelegateSimple, public LocalDynamics
    {
        Real thickness_;
        LevelSetShape *level_set_shape_;
        StdLargeVec<Vecd> &pos_, &n_, n_temp_;

      public:
        NormalPrediction(SPHBody &sph_body, Real thickness);
        virtual ~NormalPrediction(){};
        void update(size_t index_i, Real dt = 0.0);
    };

    class PredictionConvergenceCheck : public LocalDynamicsReduce<bool, ReduceAND>,
                                       public RelaxDataDelegateSimple
    {
      protected:
        const Real convergence_criterion_;
        StdLargeVec<Vecd> &n_, &n_temp_;

      public:
        PredictionConvergenceCheck(SPHBody &sph_body, Real convergence_criterion);
        virtual ~PredictionConvergenceCheck(){};

        bool reduce(size_t index_i, Real dt = 0.0);
    };

    class ConsistencyCorrection : public LocalDynamics, public RelaxDataDelegateInner
    {
      public:
        explicit ConsistencyCorrection(BaseInnerRelation &inner_relation, Real consistency_criterion);
        virtual ~ConsistencyCorrection(){};

        inline void interaction(size_t index_i, Real dt = 0.0)
        {
            mutex_modify_neighbor_.lock();
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                if (updated_indicator_[index_i] == 1)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    if (updated_indicator_[index_j] == 0)
                    {
                        updated_indicator_[index_j] = 1;
                        if (n_[index_i].dot(n_[index_j]) < consistency_criterion_)
                        {
                            if (n_[index_i].dot(-n_[index_j]) < consistency_criterion_)
                            {
                                n_[index_j] = n_[index_i];
                                updated_indicator_[index_j] = 2;
                            }
                            else
                            {
                                n_[index_j] = -n_[index_j];
                                updated_indicator_[index_j] = 1;
                            }
                        }
                    }
                }
            }
            mutex_modify_neighbor_.unlock();
        };

      protected:
        std::mutex mutex_modify_neighbor_; /**< mutex exclusion for memory conflict */
        const Real consistency_criterion_;
        StdLargeVec<int> updated_indicator_; /**> 0 not updated, 1 updated with reliable prediction, 2 updated from a reliable neighbor */
        StdLargeVec<Vecd> &n_;
    };

    class ConsistencyUpdatedCheck : public LocalDynamicsReduce<bool, ReduceAND>,
                                    public RelaxDataDelegateSimple
    {
      protected:
        StdLargeVec<int> &updated_indicator_;

      public:
        explicit ConsistencyUpdatedCheck(SPHBody &sph_body);
        virtual ~ConsistencyUpdatedCheck(){};

        bool reduce(size_t index_i, Real dt = 0.0);
    };

    class SmoothingNormal : public ParticleSmoothing<Vecd>
    {
      public:
        explicit SmoothingNormal(BaseInnerRelation &inner_relation);
        virtual ~SmoothingNormal(){};
        void update(size_t index_i, Real dt = 0.0);

      protected:
    };

    SimpleDynamics<NormalPrediction> normal_prediction_;
    ReduceDynamics<PredictionConvergenceCheck> normal_prediction_convergence_check_;
    InteractionDynamics<ConsistencyCorrection, execution::SequencedPolicy> consistency_correction_;
    ReduceDynamics<ConsistencyUpdatedCheck> consistency_updated_check_;
    InteractionWithUpdate<SmoothingNormal> smoothing_normal_;
};

/**
 * @class ShellRelaxationStepInner
 * @brief carry out particle relaxation step of particles within the shell body
 */
class ShellRelaxationStepInner : public RelaxationStepInner<PressureRelaxation>
{
  public:
    explicit ShellRelaxationStepInner(BaseInnerRelation &inner_relation, bool level_set_correction = false);
    virtual ~ShellRelaxationStepInner(){};

    SimpleDynamics<UpdateParticlePosition> update_shell_particle_position_;
    SimpleDynamics<ShellMidSurfaceBounding> mid_surface_bounding_;

    virtual void exec(Real dt = 0.0) override;
};
} // namespace relax_dynamics
} // namespace SPH
#endif // RELAX_DYNAMICS_H
