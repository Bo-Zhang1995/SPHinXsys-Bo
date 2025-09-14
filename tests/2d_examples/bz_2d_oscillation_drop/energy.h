#ifndef ENERGY_H
#define ENERGY_H

#include "complex_solid.h"
#include "elastic_dynamics.h"

namespace SPH
{
 /**
  * @class KineticEnergy
  * @brief Compute the kinematic energy
  */
	class KineticEnergy
		: public LocalDynamicsReduce<ReduceSum<Real>>
	{
	protected:
		Real *mass_;
		Vecd *vel_;

	public:
		explicit KineticEnergy(SPHBody& sph_body);
		virtual ~KineticEnergy() {};
		Real reduce(size_t index_i, Real dt = 0.0);
	};

	/**
	* @class PotentialEnergy
	* @brief Compute the Potential energy
	*/
	class PotentialEnergy
		: public LocalDynamicsReduce<ReduceSum<Real>>
	{
	private:
		SharedPtrKeeper<Gravity> gravity_ptr_keeper_;

	protected:
		Real *mass_;
		Vecd *pos_;

	public:
		explicit PotentialEnergy(SPHBody& sph_body);
		virtual ~PotentialEnergy() {};
		Real reduce(size_t index_i, Real dt = 0.0);
	};
} // namespace SPH
#endif // ENERGY_H
