#include "energy.h"

#include "base_particles.hpp"

namespace SPH
{
	//=================================================================================================//
	KineticEnergy::KineticEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<ReduceSum<Real>>(sph_body), 
		mass_(particles_->getVariableDataByName<Real>("Mass")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity"))
	{
		quantity_name_ = "KineticEnergy";
	}
	//=================================================================================================//
	Real KineticEnergy::reduce(size_t index_i, Real dt)
	{
		return 0.5 * mass_[index_i] * vel_[index_i].squaredNorm();
	}
	//=================================================================================================//
	PotentialEnergy::PotentialEnergy(SPHBody& sph_body)
		: LocalDynamicsReduce<ReduceSum<Real>>(sph_body),
		mass_(particles_->getVariableDataByName<Real>("Mass")),
        pos_(particles_->getVariableDataByName<Vecd>("Position"))
	{
		quantity_name_ = "PotentialEnergy";
	}
	//=================================================================================================//
	Real PotentialEnergy::reduce(size_t index_i, Real dt)
	{
		return 0.5 * mass_[index_i] * (pos_[index_i][0] * pos_[index_i][0] + pos_[index_i][1] * pos_[index_i][1]);
	}
	//=================================================================================================//
} // namespace SPH
//=================================================================================================//
