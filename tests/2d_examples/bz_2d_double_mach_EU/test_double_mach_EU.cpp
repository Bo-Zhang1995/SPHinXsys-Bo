/**
 * @file 	double_mach_reflection.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @author 	Bo Zhang Wang
 */
#include "common_compressible_eulerian_classes.hpp"
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                           /**< Computation domain length. */
Real DH = 1.0;                           /**< Computation domain height. */
Real particle_spacing_ref = 1.0 / 100.0; /**< Initial reference particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(DL, DH));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_one = 1.4;                         /**< initial density of one fluid. */
Real u_one = 0.0;                            /**< initial velocity of one fluid in X axis. */
Real v_one = 0.0;                            /**< initial velocity of one fluid in Y axis. */
Real p_one = 1.0;                            /**< initial pressure of one fluid. */
Real rho0_another = 8.0;                     /**< initial density of another. */
Real u_another = 8.25 * sin(3.14159 / 3.0);  /**< initial velocity of another in X axis. */
Real v_another = -8.25 * cos(3.14159 / 3.0); /**< initial velocity of another in Y axis. */
Real p_another = 140.2 / 1.2;                /**< initial pressure of another. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> CreatComputationDomain()
{
	// geometry
	std::vector<Vecd> computation_domain;
	computation_domain.push_back(Vecd(0.0, 0.0));
	computation_domain.push_back(Vecd(0.0, DH));
	computation_domain.push_back(Vecd(DL, DH));
	computation_domain.push_back(Vecd(DL, 0.0));
	computation_domain.push_back(Vecd(0.0, 0.0));
	return computation_domain;
}
class WaveBody : public ComplexShape
{
public:
	explicit WaveBody(const std::string& shape_name) : ComplexShape(shape_name)
	{
		MultiPolygon wave_block(CreatComputationDomain());
		add<MultiPolygonShape>(wave_block, "WaveBlock");
	}
};
//----------------------------------------------------------------------
//	Case-dependent initial condition.
//----------------------------------------------------------------------
class DMFInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit DMFInitialCondition(SPHBody& sph_body)
		: FluidInitialCondition(sph_body), pos_(particles_->pos_), vel_(particles_->vel_),
		rho_(particles_->rho_), p_(*particles_->getVariableByName<Real>("Pressure"))
	{
		particles_->registerVariable(mom_, "Momentum");
		particles_->registerVariable(dmom_dt_, "MomentumChangeRate");
		particles_->registerVariable(dmom_dt_prior_, "OtherMomentumChangeRate");
		particles_->registerVariable(E_, "TotalEnergy");
		particles_->registerVariable(dE_dt_, "TotalEnergyChangeRate");
		particles_->registerVariable(dE_dt_prior_, "OtherEnergyChangeRate");
		gamma_ = heat_capacity_ratio;
	};
	void update(size_t index_i, Real dt)
	{
		if (pos_[index_i][1] > tan(3.14159 / 3.0) * (pos_[index_i][0] - 1.0 / 6.0))
		{
			/** initial left wave pressure,momentum and energy profile */
			rho_[index_i] = rho0_another;
			p_[index_i] = p_another;
			Real rho_e = p_[index_i] / (gamma_ - 1.0);
			vel_[index_i][0] = u_another;
			vel_[index_i][1] = v_another;
			mom_[index_i] = rho_[index_i] * vel_[index_i];
			E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
		}
		else
		{
			rho_[index_i] = rho0_one;
			p_[index_i] = p_one;
			Real rho_e = p_[index_i] / (gamma_ - 1.0);
			vel_[index_i][0] = u_one;
			vel_[index_i][1] = v_one;
			mom_[index_i] = rho_[index_i] * vel_[index_i];
			E_[index_i] = rho_e + 0.5 * rho_[index_i] * vel_[index_i].squaredNorm();
		}
	}

protected:
	StdLargeVec<Vecd>& pos_, & vel_;
	StdLargeVec<Real>& rho_, & p_;
	StdLargeVec<Vecd> mom_, dmom_dt_, dmom_dt_prior_;
	StdLargeVec<Real> E_, dE_dt_, dE_dt_prior_;
	Real gamma_;
};


//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	// Tag for run particle relaxation for the initial body fitted distribution.
	sph_system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	sph_system.setReloadParticles(false);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody wave_block(sph_system, makeShared<WaveBody>("WaveBody"));
	wave_block.defineParticlesAndMaterial<BaseParticles, CompressibleFluid>(rho0_another, heat_capacity_ratio);
	wave_block.generateParticles<ParticleGeneratorLattice>();
	wave_block.addBodyStateForRecording<Real>("Density");
	wave_block.addBodyStateForRecording<Real>("Pressure");
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//	Note that the same relation should be defined only once.
	//----------------------------------------------------------------------
	InnerRelation wave_block_inner(wave_block);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initial condition and register variables*/
	SimpleDynamics<DMFInitialCondition> initial_condition(wave_block);
	SimpleDynamics<EulerianCompressibleTimeStepInitialization> initialize_fluid_step(wave_block);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(wave_block, 0.1);
	InteractionWithUpdate<Integration1stHalfHLLCRiemann> pressure_relaxation(wave_block_inner);
	InteractionWithUpdate<Integration2ndHalfHLLCRiemann> density_and_energy_relaxation(wave_block_inner);
	InteractionWithUpdate<KernelCorrectionMatrixInner> kernel_correction_matrix(wave_block_inner);
	InteractionDynamics<KernelGradientCorrectionInner> kernel_gradieng_update(kernel_correction_matrix);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
	RegressionTestEnsembleAverage<ReducedQuantityRecording<MaximumSpeed>>
		write_maximum_speed(io_environment, wave_block);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	initial_condition.exec();
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	write_real_body_states.writeToFile();
	kernel_correction_matrix.exec();
	kernel_gradieng_update.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 200.0; /**< End time. */
	Real D_Time = 1.0;	 /**< Time stamps for output of body states. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	write_real_body_states.writeToFile();
	write_maximum_speed.writeToFile(number_of_iterations);
	//----------------------------------------------------------------------
//	Main loop starts here.
//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < D_Time)
		{
			/** Acceleration due to viscous force. */
			initialize_fluid_step.exec();
			Real dt = get_fluid_time_step_size.exec();
			/** Dynamics including pressure relaxation. */
			integration_time += dt;
			pressure_relaxation.exec(dt);
			density_and_energy_relaxation.exec(dt);
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
					<< GlobalStaticVariables::physical_time_
					<< "	dt = " << dt << "\n";

				if (number_of_iterations % restart_output_interval == 0)
				{
					restart_io.writeToFile(number_of_iterations);
				}
			}
			number_of_iterations++;
		}

		TickCount t2 = TickCount::now();
		write_real_body_states.writeToFile();
		write_maximum_speed.writeToFile(number_of_iterations);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
