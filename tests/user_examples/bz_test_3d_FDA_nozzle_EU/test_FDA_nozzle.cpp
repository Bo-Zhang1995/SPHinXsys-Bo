/**
 * @file 	FDA_nozzle.cpp
 * @brief 	This is the 3D case for SPH Eulerian formulation.
 * @details 3D FDA nozzle example.
 * @author 	Bo Zhang
 */
#include "general_eulerian_fluid_dynamics.hpp"
#include "sphinxsys.h" 
using namespace SPH;
#define PI 3.141592653
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real throat_diameter = 0.004; // m 
Real nozzle_diameter = 0.012; // m
Real half_length = 0.1; // m
Real resolution_ref = throat_diameter / 10;
Real BW = resolution_ref * 4;
BoundingBox system_domain_bounds(
	Vec3d(-half_length, -0.5 * nozzle_diameter - BW, -0.5 * nozzle_diameter - BW),
	Vec3d(half_length, 0.5 * nozzle_diameter + BW, 0.5 * nozzle_diameter + BW));
Vecd translation(0.0, 0.0, 0.0);
Real scaling = 0.001; //transfer from m to mm.
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/FDA_nozzle_fluid.stl";
std::string full_path_to_file_wall = "./input/FDA_nozzle_wall_N10.stl";
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1056; //[kg / m3]
Real mu_f = 0.0035; //[N s / m2]
Real Q_f = 5.21e-6; 
Real U_f = Q_f * 4 / (PI * nozzle_diameter * nozzle_diameter); // Average velocity in the inlet nozzle.
Real U_max = 0.8; // From experimental data
Real c_f = U_max * 10;
//Additional information for the problem
// Re = 500 = rho_f * u_th * throat_diameter / mu_f, u_th = 0.4142992
// u_0 = u_th * (d_th/d_0)^2, u_0 = 0.04608988
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBody : public ComplexShape
{
public:
	explicit WaterBody(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
	}
};

class WallBoundary : public ComplexShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeSTL>(full_path_to_file_wall, translation, scaling);
	}
};
//----------------------------------------------------------------------
//	Application dependent initial condition
//----------------------------------------------------------------------
struct LeftInflowPressure
{
	template <class BoundaryConditionType>
	LeftInflowPressure(BoundaryConditionType& boundary_condition) {}

	Real operator()(Real& p_)
	{
		return p_;
	}
};

struct RightInflowPressure
{
	template <class BoundaryConditionType>
	RightInflowPressure(BoundaryConditionType& boundary_condition) {}

	Real operator()(Real& p_)
	{
		/*constant pressure*/
		Real pressure = 0.;
		return pressure;
	}
};
//----------------------------------------------------------------------
//	Application dependent boundary condition
//----------------------------------------------------------------------
class InletOutletBoundaryCondition
	: public LocalDynamics, public fluid_dynamics::FluidDataSimple
{
public:
	explicit InletOutletBoundaryCondition(FluidBody& fluid_body)
		: LocalDynamics(fluid_body), fluid_dynamics::FluidDataSimple(fluid_body),
		p_(*particles_->getVariableByName<Real>("Pressure")), 
		vel_(particles_->vel_), pos_(particles_->pos_){}

	void update(size_t index_i, Real dt)
	{
		Real x = pos_[index_i][0];
		Real y = pos_[index_i][1];
		Real z = pos_[index_i][2];
		Real local_radius_square = pow(y, 2) + pow(z, 2);
		if (x <= -half_length + BW)
		{
			vel_[index_i][0] = 2 * U_f * (1 - 4 * local_radius_square / nozzle_diameter / nozzle_diameter);
		}
	}
protected:
	StdLargeVec<Real>& p_;
	StdLargeVec<Vecd>& vel_, & pos_;
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class CenterLineAxialVelocityObserverParticle : public ObserverParticleGenerator
{
public:
	explicit CenterLineAxialVelocityObserverParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 51;
		Real range_of_measure = 2 * half_length;
		Real start_of_measure = -half_length;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection1VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection1VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(-0.008, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection2VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection2VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.0, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection3VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection3VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.008, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection4VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection4VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.016, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection5VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection5VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.024, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection6VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection6VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.032, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection7VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection7VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.06, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
};

class CrossSection8VelocityObserveParticle : public ObserverParticleGenerator
{
public:
	explicit CrossSection8VelocityObserveParticle(SPHBody& sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 21;
		Real range_of_measure = nozzle_diameter;
		Real start_of_measure = -0.5 * nozzle_diameter;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec3d point_coordinate(0.08, range_of_measure * (Real)i / 
				(Real)(number_of_observation_points - 1) + start_of_measure, 0.0);
			positions_.push_back(point_coordinate);
		}
	}
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
	FluidBody water_body(sph_system, makeShared<WaterBody>("WaterBody"));
	water_body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	water_body.defineAdaptationRatios(1.3, 1.0);
	water_body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? water_body.generateParticles<ParticleGeneratorReload>(io_environment, water_body.getName())
		: water_body.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	wall_boundary.defineAdaptationRatios(1.15, 1.0);
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
		: wall_boundary.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Particle and body creation of fluid observers.
	//----------------------------------------------------------------------
	ObserverBody centerlinevelocity_observer(sph_system, "CenterLineVelocityObserver");
	centerlinevelocity_observer.generateParticles<CenterLineAxialVelocityObserverParticle>();
	ObserverBody crosssection1velocity_observer(sph_system, "CrossSection1VelocityObserver");
	crosssection1velocity_observer.generateParticles<CrossSection1VelocityObserveParticle>();
	ObserverBody crosssection2velocity_observer(sph_system, "CrossSection2VelocityObserver");
	crosssection2velocity_observer.generateParticles<CrossSection2VelocityObserveParticle>();
	ObserverBody crosssection3velocity_observer(sph_system, "CrossSection3VelocityObserver");
	crosssection3velocity_observer.generateParticles<CrossSection3VelocityObserveParticle>();
	ObserverBody crosssection4velocity_observer(sph_system, "CrossSection4VelocityObserver");
	crosssection4velocity_observer.generateParticles<CrossSection4VelocityObserveParticle>();
	ObserverBody crosssection5velocity_observer(sph_system, "CrossSection5VelocityObserver");
	crosssection5velocity_observer.generateParticles<CrossSection5VelocityObserveParticle>();
	ObserverBody crosssection6velocity_observer(sph_system, "CrossSection6VelocityObserver");
	crosssection6velocity_observer.generateParticles<CrossSection6VelocityObserveParticle>();
	ObserverBody crosssection7velocity_observer(sph_system, "CrossSection7VelocityObserver");
	crosssection7velocity_observer.generateParticles<CrossSection7VelocityObserveParticle>();
	ObserverBody crosssection8velocity_observer(sph_system, "CrossSection8VelocityObserver");
	crosssection8velocity_observer.generateParticles<CrossSection8VelocityObserveParticle>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation water_block_inner(water_body);
	InnerRelation wall_boundary_inner(wall_boundary);
	ContactRelation water_wall_contact(water_body, { &wall_boundary });
	ComplexRelation water_block_complex(water_body, { &wall_boundary });
	ComplexRelation water_block_complex_corrected(water_body, { &wall_boundary });
	ComplexRelation wall_boundary_block_complex(wall_boundary, { &water_body });//only used for relaxation
	//----------------------------------------------------------------------
	//  Observer relation
	//---------------------------------------------------------------------- 
	ContactRelation centerlinevelocity_observer_contact(centerlinevelocity_observer, { &water_body });
	ContactRelation crosssection1velocity_observer_contact(crosssection1velocity_observer, { &water_body });
	ContactRelation crosssection2velocity_observer_contact(crosssection2velocity_observer, { &water_body });
	ContactRelation crosssection3velocity_observer_contact(crosssection3velocity_observer, { &water_body });
	ContactRelation crosssection4velocity_observer_contact(crosssection4velocity_observer, { &water_body });
	ContactRelation crosssection5velocity_observer_contact(crosssection5velocity_observer, { &water_body });
	ContactRelation crosssection6velocity_observer_contact(crosssection6velocity_observer, { &water_body });
	ContactRelation crosssection7velocity_observer_contact(crosssection7velocity_observer, { &water_body });
	ContactRelation crosssection8velocity_observer_contact(crosssection8velocity_observer, { &water_body });
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_body);
		SimpleDynamics<RandomizeParticlePosition> random_boundary_particles(wall_boundary);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_water_body_to_file(io_environment, { &water_body });
		BodyStatesRecordingToVtp write_boundary_body_to_file(io_environment, { &wall_boundary });
		/** Write the particle reload files. */
		ReloadParticleIO write_water_particle_reload(io_environment, water_body);
		ReloadParticleIO write_boundary_particle_reload(io_environment, wall_boundary);
		/* Relaxation methods including pressure-based and KGC-based relaxation.*/
		InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex_water(water_block_complex);
		InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_complex_column(wall_boundary_block_complex);
		relax_dynamics::RelaxationStepComplexImplicit relaxation_water_inner(water_block_complex);
		relax_dynamics::RelaxationStepInnerImplicit relaxation_boundary_inner(wall_boundary_inner, true);
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_water_block_kinetic_energy(water_block_inner);
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_wall_boundary_kinetic_energy(wall_boundary_inner);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_water_block_average_kinetic_energy(water_body, "ParticleKineticEnergy");
		ReduceDynamics<QuantityMaximum<Real>> calculate_water_block_maximum_kinetic_energy(water_body, "ParticleKineticEnergy");
		water_body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
		//----------------------------------------------------------------------
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_water_body_particles.exec(0.25);
		random_boundary_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		relaxation_water_inner.SurfaceBounding().exec();
		relaxation_boundary_inner.SurfaceBounding().exec();
		write_boundary_body_to_file.writeToFile(0);
		write_water_body_to_file.writeToFile(0);
		Real water_block_average_energy = 100.0;
		Real water_block_maximum_energy = 100.0;
		Real last_water_block_maximum_energy = 100.0;
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.
		while (ite < 20000)
		{
			relaxation_boundary_inner.exec();
			relaxation_water_inner.exec();

			if (ite % 2000 == 0)
			{
				update_water_block_kinetic_energy.exec();
				water_block_average_energy = calculate_water_block_average_kinetic_energy.exec();
				water_block_maximum_energy = calculate_water_block_maximum_kinetic_energy.exec();

				std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
				std::cout << "Body: " << " Maximum: " << water_block_maximum_energy << " Average: " 
					<< water_block_average_energy << std::endl;

				write_water_body_to_file.writeToFile(ite);
				write_boundary_body_to_file.writeToFile(ite);
			}
			ite++;
		}

		update_wall_boundary_kinetic_energy.exec();
		std::cout << "Residual: " << " maximum: " << calculate_water_block_maximum_kinetic_energy.exec()
			<< " average: " << calculate_water_block_average_kinetic_energy.exec() << std::endl;
		std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;
		write_water_body_to_file.writeToFile(ite);
		write_boundary_body_to_file.writeToFile(ite);
		write_water_particle_reload.writeToFile(0);
		write_boundary_particle_reload.writeToFile(0);

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
	SimpleDynamics<NozzleInletInitialCondition> nozzle_inlet_initial_condition(wall_boundary);
	SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
	SimpleDynamics<InletOutletBoundaryCondition> nozzle_inlet_outlet_condition(water_body);
	/** Intialize a time step */
	SimpleDynamics<TimeStepInitialization> time_step_initialization(water_body);
	/** Initialize the correction matrix and correctecd configuration */
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_matrix(water_block_complex);
	InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_matrix_corrected(water_block_complex_corrected);
	InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(kernel_correction_matrix_corrected);
	/** Time step size with considering sound wave speed */
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_body);
	InteractionDynamics<fluid_dynamics::DistanceFromWall> fluid_distance_to_wall(water_block_complex); //fluid particle distance to wall/
	InteractionDynamics<fluid_dynamics::BoundingFromWall> fluid_bounding_distance(water_block_complex);
	InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex_corrected);
	/** Pressure relaxation algorithm by using verlet time stepping */
	InteractionWithUpdate<fluid_dynamics::ICEIntegration1stHalfHLLERiemannWithWall> pressure_relaxation(water_block_complex);
	InteractionWithUpdate<fluid_dynamics::ICEIntegration2ndHalfHLLERiemannWithWall> density_and_energy_relaxation(water_block_complex);
	water_body.addBodyStateForRecording<Real>("Density");
	water_body.addBodyStateForRecording<Vecd>("Position");
	water_body.addBodyStateForRecording<Real>("Pressure");
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");
	wall_boundary.addBodyStateForRecording<Vecd>("Position");
	//----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //---------------------------------------------------------------------- 
	/** Output the body states */
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	/** Output the body states for restart simulation */
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	/** Output the observed velocity */
	ObservedQuantityRecording<Vecd> write_centerline_velocity("Velocity", io_environment, centerlinevelocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection1_velocity("Velocity", io_environment, crosssection1velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection2_velocity("Velocity", io_environment, crosssection2velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection3_velocity("Velocity", io_environment, crosssection3velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection4_velocity("Velocity", io_environment, crosssection4velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection5_velocity("Velocity", io_environment, crosssection5velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection6_velocity("Velocity", io_environment, crosssection6velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection7_velocity("Velocity", io_environment, crosssection7velocity_observer_contact);
	ObservedQuantityRecording<Vecd> write_crosssection8_velocity("Velocity", io_environment, crosssection8velocity_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	nozzle_inlet_initial_condition.exec();
	wall_boundary_normal_direction.exec();
	fluid_distance_to_wall.exec();
	fluid_bounding_distance.exec();
	kernel_correction_matrix.exec();
	kernel_correction_matrix_corrected.exec();
	kernel_gradient_update.exec();
	body_states_recording.writeToFile();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = 0;
	int screen_output_interval = 100;
	int restart_output_interval = screen_output_interval * 10;
	Real End_Time = 5; /**< End time. */
	Real D_Time = End_Time / 200;	 /**< Time stamps for output of body states. */
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
			nozzle_inlet_outlet_condition.exec();
			density_and_energy_relaxation.exec(dt);
			GlobalStaticVariables::physical_time_ += dt;

			if (number_of_iterations % screen_output_interval == 0)
			{
				std::cout << std::fixed << std::setprecision(9) << "N=" 
					<< number_of_iterations << "	Time = "
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
		body_states_recording.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
		
		write_centerline_velocity.writeToFile(number_of_iterations);
		write_crosssection1_velocity.writeToFile(number_of_iterations);
		write_crosssection2_velocity.writeToFile(number_of_iterations);
		write_crosssection3_velocity.writeToFile(number_of_iterations);
		write_crosssection4_velocity.writeToFile(number_of_iterations);
		write_crosssection5_velocity.writeToFile(number_of_iterations);
		write_crosssection6_velocity.writeToFile(number_of_iterations);
		write_crosssection7_velocity.writeToFile(number_of_iterations);
		write_crosssection8_velocity.writeToFile(number_of_iterations);
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}