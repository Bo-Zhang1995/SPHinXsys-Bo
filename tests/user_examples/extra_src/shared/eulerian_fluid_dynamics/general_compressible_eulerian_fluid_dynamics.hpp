/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                         *
 * ------------------------------------------------------------------------*/
#pragma once

#include "general_compressible_eulerian_fluid_dynamics.h"

namespace SPH
{
//=================================================================================================//
template <class RiemannSolverType>
CEIntegration1stHalf<RiemannSolverType>::CEIntegration1stHalf(BaseInnerRelation& inner_relation)
    : BaseIntegrationInGeneralCE(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_) {};
//=================================================================================================//
template <class RiemannSolverType>
void CEIntegration1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
    Vecd momentum_change_rate = dmom_dt_prior_[index_i];
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd& e_ij = inner_neighborhood.e_ij_[n];

        CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);

        Vec2d wave_speeds = riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
        Real s_l = wave_speeds[0];
        Real s_r = wave_speeds[1];

        Matd flux_l = (rho_[index_i] * vel_[index_i] * vel_[index_i].transpose() + p_[index_i] * Matd::Identity()) * B_[index_j];
        Matd flux_r = (rho_[index_j] * vel_[index_j] * vel_[index_j].transpose() + p_[index_j] * Matd::Identity()) * B_[index_i];

        if (s_l < 0 && s_r > 0)
        {
            DissipationStateCE dissipation_state = riemann_solver_.getDissipationState(state_i, state_j, e_ij);
            momentum_change_rate -= 2 * (0.5 * (flux_l + flux_r) + dissipation_state.momentum_dissipation_) * e_ij * dW_ijV_j;
        }
        else if (s_l > 0)
        {
            momentum_change_rate -= 2 * flux_l * e_ij * dW_ijV_j;
        }
        else if (s_r < 0)
        {
            momentum_change_rate -= 2 * flux_r * e_ij * dW_ijV_j;
        }
    }
    dmom_dt_[index_i] = momentum_change_rate;
}
//=================================================================================================//
template <class RiemannSolverType>
void CEIntegration1stHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    mom_[index_i] += dmom_dt_[index_i] * dt;
    vel_[index_i] = mom_[index_i] / rho_[index_i];
}
//=================================================================================================//
template <class RiemannSolverType>
CEIntegration2ndHalf<RiemannSolverType>::CEIntegration2ndHalf(BaseInnerRelation& inner_relation)
    : BaseIntegrationInGeneralCE(inner_relation), riemann_solver_(compressible_fluid_, compressible_fluid_) {};
//=================================================================================================//
template <class RiemannSolverType>
void CEIntegration2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
{
    CompressibleFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], E_[index_i]);
    Real density_change_rate = 0.0;
    Real energy_change_rate = dE_dt_prior_[index_i];
    Neighborhood& inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real dW_ijV_j = inner_neighborhood.dW_ijV_j_[n];
        Vecd& e_ij = inner_neighborhood.e_ij_[n];

        CompressibleFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], E_[index_j]);

        Vec2d wave_speeds = riemann_solver_.getBoundingWaveSpeeds(state_i, state_j, e_ij);
        Real s_l = wave_speeds[0];
        Real s_r = wave_speeds[1];

        Vecd flux_l = rho_[index_i] * vel_[index_i] * e_ij.transpose() * B_[index_j] * e_ij;
        Vecd flux_r = rho_[index_j] * vel_[index_j] * e_ij.transpose() * B_[index_i] * e_ij;

        Vecd flux_energy_l = (E_[index_i] * vel_[index_i] + p_[index_i] * vel_[index_i]) * e_ij.transpose() * B_[index_j] * e_ij;
        Vecd flux_energy_r = (E_[index_j] * vel_[index_j] + p_[index_j] * vel_[index_j]) * e_ij.transpose() * B_[index_i] * e_ij;

        if (s_l < 0 && s_r > 0)
        {
            DissipationStateCE dissipation_state = riemann_solver_.getDissipationState(state_i, state_j, e_ij);
            density_change_rate -= 2 * ((0.5 * (flux_l + flux_r) + dissipation_state.density_dissipation_).dot(e_ij)) * dW_ijV_j;
            energy_change_rate -= 2 * ((0.5 * (flux_energy_l + flux_energy_r) + dissipation_state.energy_dissipation_).dot(e_ij)) * dW_ijV_j;
        }
        else if (s_l > 0)
        {
            density_change_rate -= 2 * flux_l.dot(e_ij) * dW_ijV_j;
            energy_change_rate -= 2 * flux_energy_l.dot(e_ij) * dW_ijV_j;
        }
        else if (s_r < 0)
        {
            density_change_rate -= 2 * flux_r.dot(e_ij) * dW_ijV_j;
            energy_change_rate -= 2 * flux_energy_r.dot(e_ij) * dW_ijV_j;
        }
    }
    drho_dt_[index_i] = density_change_rate;
    dE_dt_[index_i] = energy_change_rate;
};
//=================================================================================================//
template <class RiemannSolverType>
void CEIntegration2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
{
    E_[index_i] += dE_dt_[index_i] * dt;
    rho_[index_i] += drho_dt_[index_i] * dt;
    Real rho_e = E_[index_i] - 0.5 * mom_[index_i].squaredNorm() / rho_[index_i];
    p_[index_i] = compressible_fluid_.getPressure(rho_[index_i], rho_e);
}
//=================================================================================================//
} // namespace SPH
//=================================================================================================//