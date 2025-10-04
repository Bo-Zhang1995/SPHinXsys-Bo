#pragma once

#include "fluid_dynamics_inner_correction.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration1stHalfConsistency<RiemannSolverType>::BaseIntegration1stHalfConsistency(BaseInnerRelation& inner_relation)
    : BaseIntegration1stHalf<RiemannSolverType>(inner_relation),
    B_(*this->particles_->template registerSharedVariable<Matd>("KernelCorrectionMatrix", Matd::Identity())) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalfConsistency<RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    BaseIntegration1stHalf<RiemannSolverType>::initialization(index_i, dt);
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration1stHalfConsistency<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Vecd acceleration = Vecd::Zero();
    Real rho_dissipation(0);
    const Neighborhood& inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        const Vecd& e_ij = inner_neighborhood.e_ij_[n];

        acceleration -= (this->p_[index_i] * this->B_[index_j] + this->p_[index_j] * this->B_[index_i]) * dW_ijV_j * e_ij;
        rho_dissipation += this->riemann_solver_.DissipativeUJump(this->p_[index_i] - this->p_[index_j]) * dW_ijV_j;
    }
    this->acc_[index_i] += acceleration  / this->rho_[index_i];
    this->drho_dt_[index_i] = rho_dissipation * this->rho_[index_i];
    this->vel_div_[index_i] = abs(rho_dissipation);
}
//=================================================================================================//
template <class RiemannSolverType>
BaseIntegration2ndHalfCorrect<RiemannSolverType>::BaseIntegration2ndHalfCorrect(BaseInnerRelation &inner_relation)
    : BaseIntegration2ndHalf<RiemannSolverType>(inner_relation),
      B_(*this->particles_->template registerSharedVariable<Matd>("KernelCorrectionMatrix", Matd::Identity())) {}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfCorrect<RiemannSolverType>::initialization(size_t index_i, Real dt)
{
    BaseIntegration2ndHalf<RiemannSolverType>::initialization(index_i, dt);
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfCorrect<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    Real density_change_rate(0);
    Vecd p_dissipation = Vecd::Zero();
    const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        const Vecd &e_ij = inner_neighborhood.e_ij_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];

        Real u_jump = (this->vel_[index_i] - this->vel_[index_j]).dot(e_ij);
        density_change_rate += (this->vel_[index_i] - this->vel_[index_j]).dot(0.5 * (B_[index_i] + B_[index_j]) * e_ij) * dW_ijV_j;
        p_dissipation += this->riemann_solver_.DissipativePJump(u_jump) * dW_ijV_j * e_ij;
    }
    this->drho_dt_[index_i] += density_change_rate * this->rho_[index_i];
    this->vel_div_[index_i] += abs(density_change_rate);
    this->acc_[index_i] = p_dissipation / this->rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
void BaseIntegration2ndHalfCorrect<RiemannSolverType>::update(size_t index_i, Real dt)
{
    BaseIntegration2ndHalf<RiemannSolverType>::update(index_i, dt);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
