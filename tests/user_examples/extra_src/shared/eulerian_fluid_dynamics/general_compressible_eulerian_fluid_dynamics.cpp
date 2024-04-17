
#include "general_compressible_eulerian_fluid_dynamics.h"

namespace SPH
{
//=================================================================================================//
Vec2d NoRiemannSolverCE::getBoundingWaveSpeeds(const CompressibleFluidState& state_i, 
    const CompressibleFluidState& state_j, const Vecd e_ij)
{
    Real s_l = -1;
    Real s_r = 1;
    return Vec2d(s_l, s_r);
}
//=================================================================================================//
DissipationStateCE NoRiemannSolverCE::getDissipationState(const CompressibleFluidState& state_i, 
    const CompressibleFluidState& state_j, const Vecd e_ij)
{
    return DissipationStateCE(); //default return zero.
}
//=================================================================================================//
Vec2d HLLERiemannSolverCE::getBoundingWaveSpeeds(const CompressibleFluidState& state_i, const
    CompressibleFluidState& state_j, const Vecd e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real R_lf = state_j.rho_ / state_i.rho_;
    Real u_tlide = (ul + ur * R_lf) / (1.0 + R_lf);
    Real v_tlide = ((state_i.vel_ - ul * (-e_ij)).norm() + (state_j.vel_ - ur * (-e_ij)).norm() * R_lf) / (1.0 + R_lf);
    Real hl = (state_i.E_ + state_i.p_) / state_i.rho_;
    Real hr = (state_j.E_ + state_j.p_) / state_j.rho_;
    Real h_tlide = (hl + hr * R_lf) / (1.0 + R_lf);
    Real sound_tlide = sqrt((1.4 - 1.0) * (h_tlide - 0.5 * (u_tlide * u_tlide + v_tlide * v_tlide)));
    Real s_l = SMIN(ul - compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_), u_tlide - sound_tlide);
    Real s_r = SMAX(ur + compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_), u_tlide + sound_tlide);
    return Vec2d(s_l, s_r);
};
//=================================================================================================//
DissipationStateCE HLLERiemannSolverCE::getDissipationState(const CompressibleFluidState& state_i, 
    const CompressibleFluidState& state_j, const Vecd e_ij)
{
    Real ul = -e_ij.dot(state_i.vel_);
    Real ur = -e_ij.dot(state_j.vel_);
    Real cl = this->compressible_fluid_i_.getSoundSpeed(state_i.p_, state_i.rho_);
    Real cr = this->compressible_fluid_j_.getSoundSpeed(state_j.p_, state_j.rho_);
    Real clr = (cl * state_i.rho_ + cr * state_j.rho_) / (state_i.rho_ + state_j.rho_);

    Real s_l = getBoundingWaveSpeeds(state_i, state_j, e_ij)[0];
    Real s_r = getBoundingWaveSpeeds(state_i, state_j, e_ij)[1];

    Matd flux_lm = (state_i.rho_ * state_i.vel_ * state_i.vel_.transpose() + state_i.p_ * Matd::Identity());
    Matd flux_rm = (state_j.rho_ * state_j.vel_ * state_j.vel_.transpose() + state_j.p_ * Matd::Identity());
    Vecd state_lm = state_i.rho_ * state_i.vel_;
    Vecd state_rm = state_j.rho_ * state_j.vel_;

    Vecd flux_ld = (state_i.rho_ * state_i.vel_);
    Vecd flux_rd = (state_j.rho_ * state_j.vel_);
    Real state_ld = state_i.rho_;
    Real state_rd = state_j.rho_;

    Vecd flux_le = (state_i.E_ * state_i.vel_ + state_i.p_ * state_i.vel_);
    Vecd flux_re = (state_j.E_ * state_j.vel_ + state_j.p_ * state_j.vel_);
    Real state_le = state_i.E_;
    Real state_re = state_j.E_;
    
    //Matd momentum_dissipation = SMIN<Real>(20.0 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) * 
    //    (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_lm - flux_rm) + s_r * s_l * (state_rm - state_lm) * (-e_ij).transpose() / (s_r - s_l));
    //Vecd density_dissipation = SMIN<Real>(20.0 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) * 
    //    (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_ld - flux_rd) + s_r * s_l * (state_rd - state_ld) / (s_r - s_l) * (-e_ij));
    //Vecd energy_dissipation = SMIN<Real>(20.0 * SMAX<Real>((ul - ur) / clr, Real(0)), Real(1)) *
    //    (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_le - flux_re) + s_r * s_l * (state_re - state_le) / (s_r - s_l) * (-e_ij));

    Matd momentum_dissipation =
        (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_lm - flux_rm) + s_r * s_l * (state_rm - state_lm) * (-e_ij).transpose() / (s_r - s_l));
    Vecd density_dissipation = 
        (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_ld - flux_rd) + s_r * s_l * (state_rd - state_ld) / (s_r - s_l) * (-e_ij));
    Vecd energy_dissipation = 
        (0.5 * (s_r + s_l) / (s_r - s_l) * (flux_le - flux_re) + s_r * s_l * (state_re - state_le) / (s_r - s_l) * (-e_ij));

    return DissipationStateCE(momentum_dissipation, density_dissipation, energy_dissipation);
};
//=================================================================================================//
EulerianCEViscousAccelerationInner::EulerianCEViscousAccelerationInner(BaseInnerRelation& inner_relation)
    : ViscousAccelerationInner(inner_relation),
    dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
    dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")) {};
//=================================================================================================//
void EulerianCEViscousAccelerationInner::interaction(size_t index_i, Real dt)
{
    Real rho_i = rho_[index_i];
    const Vecd& vel_i = vel_[index_i];

    Vecd acceleration = Vecd::Zero();
    Vecd vel_derivative = Vecd::Zero();
    const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];

        // viscous force
        vel_derivative = (vel_i - vel_[index_j]) / (inner_neighborhood.r_ij_[n] + 0.01 * smoothing_length_);
        acceleration += 2.0 * mu_ * vel_derivative * inner_neighborhood.dW_ijV_j_[n] / rho_i;
    }
    dmom_dt_prior_[index_i] += rho_[index_i] * acceleration;
    dE_dt_prior_[index_i] += rho_[index_i] * acceleration.dot(vel_[index_i]);
}
//=================================================================================================//
EulerianCETimeStepInitialization::EulerianCETimeStepInitialization(SPHBody& sph_body, SharedPtr<Gravity> gravity_ptr)
    : TimeStepInitialization(sph_body, gravity_ptr), rho_(particles_->rho_), 
    pos_(particles_->pos_), vel_(particles_->vel_),
    dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")),
    dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")) {};
//=================================================================================================//
void EulerianCETimeStepInitialization::update(size_t index_i, Real dt)
{
    dmom_dt_prior_[index_i] = rho_[index_i] * gravity_->InducedAcceleration(pos_[index_i]);
    dE_dt_prior_[index_i] = rho_[index_i] * (gravity_->InducedAcceleration(pos_[index_i])).dot(vel_[index_i]);
}
//=================================================================================================//
EulerianCEAcousticTimeStepSize::EulerianCEAcousticTimeStepSize(SPHBody& sph_body, Real acousticCFL)
    : AcousticTimeStepSize(sph_body), rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure")), 
    vel_(particles_->vel_),
    smoothing_length_(sph_body.sph_adaptation_->ReferenceSmoothingLength()),
    compressible_fluid_(CompressibleFluid(1.0, 1.4)), acousticCFL_(acousticCFL) {};
//=================================================================================================//
Real EulerianCEAcousticTimeStepSize::reduce(size_t index_i, Real dt)
{
    return compressible_fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real EulerianCEAcousticTimeStepSize::outputResult(Real reduced_value)
{
    return acousticCFL_ / Dimensions * smoothing_length_ / (reduced_value + TinyReal);
}
//=================================================================================================//
BaseIntegrationInGeneralCE::BaseIntegrationInGeneralCE(BaseInnerRelation& inner_relation)
    : BaseIntegration(inner_relation), compressible_fluid_(CompressibleFluid(1.0, 1.4)),
    Vol_(particles_->Vol_), E_(*particles_->getVariableByName<Real>("TotalEnergy")),
    dE_dt_(*particles_->getVariableByName<Real>("TotalEnergyChangeRate")),
    dE_dt_prior_(*particles_->getVariableByName<Real>("OtherEnergyChangeRate")),
    mom_(*particles_->getVariableByName<Vecd>("Momentum")),
    dmom_dt_(*particles_->getVariableByName<Vecd>("MomentumChangeRate")),
    dmom_dt_prior_(*particles_->getVariableByName<Vecd>("OtherMomentumChangeRate")) {};
//=================================================================================================//
}// namespace SPH
