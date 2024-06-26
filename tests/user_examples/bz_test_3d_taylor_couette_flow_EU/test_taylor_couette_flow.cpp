/**
 * @file 	taylor_couette_flow.cpp
 * @brief 	This is the one of the basic test cases for SPH Eulerian formulation.
 * @details 3D eulerian_taylor_couette_flow example.
 * @author 	Bo Zhang
 */
#include "general_eulerian_fluid_dynamics.hpp" //eulerian classes for fluid.
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Setting for the second geometry.
//	To use this, please commenting the setting for the first geometry.
//----------------------------------------------------------------------
std::string full_path_to_file_fluid = "./input/case.stl";
std::string full_path_to_file_wall = "./input/wall.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real R_in = 1.5;   /* Inner radius. */
Real R_out = 2.5;  /* Outer radius. */
Real L = 4.2;
Real d = R_out - R_in; /* the radius difference. */
Real resolution_ref = d / 40;
Real BW = resolution_ref * 5;
Real distance = 0.2;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-R_out - distance, -R_out - distance, -0.5 * L - distance), Vec3d(R_out + distance, R_out + distance, 0.5 * L + distance));
int resolution(10);
Vec3d translation(0.0, 0.0, 0.0);
Real scaling = 1;
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;					        /**< Reference density of fluid. */
Real U_f = 0.1;                             /**< Angular velocity. */
Real U_max = U_f *  R_in;                   /**< Reference maximum velocity. */
Real c_f = 10 * U_max;                      /**< Reference sound speed. */
Real Re = 150;                              /**< Reynolds number. */
Real mu_f = rho0_f * U_f * R_in / Re;       /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class FluidColumn : public ComplexShape
{
public:
	explicit FluidColumn(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_file_fluid, translation, scaling);
	}
};

class OuterColumn : public ComplexShape
{
public:
	explicit OuterColumn(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_file_wall, translation, scaling);
	}
};

//----------------------------------------------------------------------
//	Application dependent initial condition
//----------------------------------------------------------------------
class RotatingColumnInitialCondition
	: public LocalDynamics, public SolidDataSimple
{
public:
	explicit RotatingColumnInitialCondition(SolidBody& solid_body)
		: LocalDynamics(solid_body), SolidDataSimple(solid_body),
		vel_(particles_->vel_), pos_(particles_->pos_) {};

	void update(size_t index_i, Real dt)
	{
		Real x = pos_[index_i][0];
		Real y = pos_[index_i][1];
		Real z = pos_[index_i][2];
		Real angular_velocity = U_f;
		Real local_radius = sqrt(pow(x, 2) + pow(y, 2));
		Real angular = atan2(y, x);

		if ((z >= -2.5) && (z <= 2.5) && (sqrt(x * x + y * y) <= 1.5))
		{
			vel_[index_i][0] = angular_velocity * local_radius * sin(angular);
			vel_[index_i][1] = -angular_velocity * local_radius * cos(angular);
		}
	}
protected:
	StdLargeVec<Vecd>& vel_, & pos_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	// Tag for run particle relaxation for tHe initial body fitted distribution.
	sph_system.setRunParticleRelaxation(false);
	// Tag for computation start with relaxed body fitted particles distribution.
	sph_system.setReloadParticles(true);
	/** Set the starting time. */
	GlobalStaticVariables::physical_time_ = 0.0;
	IOEnvironment io_environment(sph_system);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody water_body(sph_system, makeShared<FluidColumn>("WaterBody"));
	water_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	water_body.defineAdaptationRatios(0.95, 1.0);
	water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	//water_body.generateParticles<ParticleGeneratorLattice>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
		: water_body.generateParticles<ParticleGeneratorLattice>();
	water_body.addBodyStateForRecording<Real>("Density");
	water_body.addBodyStateForRecording<Vecd>("Position");
	/**
	 * @brief 	Particle and body creation of wall boundary.
	 */
	SolidBody outer_column(sph_system, makeShared<OuterColumn>("OuterColumn"));
	outer_column.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	outer_column.defineAdaptationRatios(1.15, 1.0);
	outer_column.defineParticlesAndMaterial<SolidParticles, Solid>();
	//outer_column.generateParticles<ParticleGeneratorLattice>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? outer_column.generateParticles<ParticleGeneratorReload>(io_environment, outer_column.getName())
		: outer_column.generateParticles<ParticleGeneratorLattice>();
	outer_column.addBodyStateForRecording<Vecd>("NormalDirection");
	outer_column.addBodyStateForRecording<Vecd>("Position");
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_body);
	ComplexRelation water_block_complex(water_body, { &outer_column });
	ComplexRelation water_block_complex_corrected(water_body, { &outer_column });
	ComplexRelation column_block_complex(outer_column, { &water_body });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		/** body topology only for particle relaxation */
		InnerRelation water_block_inner(water_body);
		InnerRelation boundary_inner(outer_column);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_body);
		SimpleDynamics<RandomizeParticlePosition> random_boundary_particles(outer_column);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToPlt write_water_body_to_vtp(io_environment, { &water_body });
		BodyStatesRecordingToPlt write_boundary_body_to_vtp(io_environment, { &outer_column });
		/** Write the particle reload files. */
		ReloadParticleIO write_water_particle_reload(io_environment, water_body);
		ReloadParticleIO write_boundary_particle_reload(io_environment, outer_column);

		/* Relaxation method: including based on the 0th and 1st order consistency. */
		InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex_water(water_block_complex);
		InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex_column(column_block_complex);
		relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_water_inner(water_block_inner, true);
		relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_boundary_inner(boundary_inner, true);
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_water_block_kinetic_energy(water_block_inner);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_water_block_average_kinetic_energy(water_body, "ParticleKineticEnergy");
		ReduceDynamics<QuantityMaximum<Real>> calculate_water_block_maximum_kinetic_energy(water_body, "ParticleKineticEnergy");
		//water_body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_water_body_particles.exec(0.25);
		random_boundary_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		relaxation_water_inner.SurfaceBounding().exec();
		relaxation_boundary_inner.SurfaceBounding().exec();
		write_water_body_to_vtp.writeToFile(0);
		Real water_block_average_energy = 100.0;
		Real water_block_maximum_energy = 100.0;
		Real last_water_block_maximum_energy = 100.0;
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.
		GlobalStaticVariables::physical_time_ = ite;
		while (ite < 10000)
		{
			kernel_correction_complex_column.exec();
			relaxation_boundary_inner.exec();
			kernel_correction_complex_water.exec();
			relaxation_water_inner.exec();

			ite += 1;
			if (ite % 1000 == 0)
			{
				update_water_block_kinetic_energy.exec();
				water_block_average_energy = calculate_water_block_average_kinetic_energy.exec();
				water_block_maximum_energy = calculate_water_block_maximum_kinetic_energy.exec();
				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				std::cout << "Body: " << " Maximum: " << water_block_maximum_energy << " Average: " << water_block_average_energy << std::endl;

				write_water_body_to_vtp.writeToFile(ite);
				write_boundary_body_to_vtp.writeToFile(ite);
			}
		}

		std::cout << "Residual: " << " maximum: " << calculate_water_block_maximum_kinetic_energy.exec() << " average: " 
			<< calculate_water_block_average_kinetic_energy.exec() << std::endl;
		std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

		ite++;
		write_water_body_to_vtp.writeToFile(ite);
		write_boundary_body_to_vtp.writeToFile(ite);
		write_boundary_particle_reload.writeToFile(0);
		write_water_particle_reload.writeToFile(0);

		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	/** Initial condition with momentum and energy field */
	SimpleDynamics<RotatingColumnInitialCondition> solid_initial_condition(outer_column);
	SimpleDynamics<NormalDirectionFromBodyShape> boundary_normal_direction_outer(outer_column);
	/** Initialize a time step. */
	SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_maxtrix(water_block_complex);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_maxtrix_corrected(water_block_complex_corrected);
	InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(kernel_correction_maxtrix_corrected);
	/** Time step size with considering sound wave speed. */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_corrected);
	/** Pressure relaxation algorithm by using verlet time stepping. */
	InteractionWithUpdate<fluid_dynamics::ICEIntegration1stHalfHLLERiemannWithWall> pressure_relaxation(water_block_complex);
	InteractionWithUpdate<fluid_dynamics::ICEIntegration2ndHalfHLLERiemannWithWall> density_and_energy_relaxation(water_block_complex);
	/** Computing vorticity in the flow. */
	InteractionDynamics<fluid_dynamics::VorticityWithWall> compute_vorticity_xyz(water_block_complex_corrected);
	InteractionWithUpdate<fluid_dynamics::AngleVorticityWithWall> compute_vorticity(water_block_complex_corrected);
	water_body.addBodyStateForRecording<Real>("ThetaVorticity");
	water_body.addBodyStateForRecording<Mat3d>("VelocityGradient");
	water_body.addBodyStateForRecording<AngularVecd>("Vorticity");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	/** Output the body states. */
	BodyStatesRecordingToPlt body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the body states for restart simulation. */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	solid_initial_condition.exec();
	boundary_normal_direction_outer.exec();
	kernel_correction_maxtrix.exec();
	kernel_correction_maxtrix_corrected.exec();
	kernel_gradient_update.exec();
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 1000.0; /**< End time. */
	Real D_Time = 10.0;	 /**< Time stamps for output of body states. */
	/** statistics for computing CPU time. */
	TickCount t1 = TickCount::now();
    TimeInterval interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	/** Output the start states of bodies. */
	body_states_recording.writeToFile();
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
			time_step_initialization.exec();
			Real dt = get_fluid_time_step_size.exec();
			viscous_acceleration.exec();
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
		compute_vorticity.exec();
		compute_vorticity_xyz.exec();
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
        interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
