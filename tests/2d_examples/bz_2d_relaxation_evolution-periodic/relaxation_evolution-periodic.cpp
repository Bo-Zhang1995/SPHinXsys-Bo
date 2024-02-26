/**
 * @file 	relaxation_evolution_periodic.cpp
 * @brief   This is the first case with periodical boundary condition
 *          by testing the relaxation and consistency.
 * @author  Bo Zhang, Xiangyu Hu.
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;
Real LH = 1.0;
Real resolution_ref = LH / 40.0;
Real BW = resolution_ref * 4;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(LL + BW, LH + BW));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d water_block_translation = water_block_halfsize;

class Body : public ComplexShape
{
public:
	explicit Body(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform(water_block_halfsize), water_block_translation);
	}
};

int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true);
	sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody body(sph_system, makeShared<Body>("RelaxationBody"));
	body.defineAdaptationRatios(1.3, 1.0);
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1, 1, 1);
	body.addBodyStateForRecording<Vecd>("Position");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
		: body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map.
		//----------------------------------------------------------------------
		InnerRelation body_inner(body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Plt file. */
		BodyStatesRecordingToPlt write_body_to_plt(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });
		/* Relaxation method : including based on the 0th and 1st order consistency. */
		InteractionWithUpdate<KernelCorrectionMatrixInner> correction_matrix(body_inner);
		body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
		//relax_dynamics::RelaxationStepInner relaxation_inner(body_inner, false);
		relax_dynamics::RelaxationStepInnerImplicit relaxation_inner(body_inner, false);
		//relax_dynamics::RelaxationStepInner<CorrectionMatrixRelaxation> relaxation_inner(body_inner, false);
		//relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner(body_inner, false);
		/* Periodic boundary condition. */
		PeriodicConditionUsingCellLinkedList periodic_condition_x(body, body.getBodyShapeBounds(), xAxis);
		PeriodicConditionUsingCellLinkedList periodic_condition_y(body, body.getBodyShapeBounds(), yAxis);
		/** Update the kinetic energy for statistical. */
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(body_inner);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_kinetic_energy(body, "ParticleKineticEnergy");
		ReduceDynamics<QuantityMaximum<Real>> calculate_body_maximum_kinetic_energy(body, "ParticleKineticEnergy");
		std::string filefullpath_error_analysis = io_environment.output_folder_ + "/" + "error_analysis.dat";
		std::ofstream out_file_error_analysis(filefullpath_error_analysis.c_str(), std::ios::app);
		//----------------------------------------------------------------------  
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		periodic_condition_x.update_cell_linked_list_.exec();
		periodic_condition_y.update_cell_linked_list_.exec();
		sph_system.initializeSystemConfigurations();
		write_body_to_plt.writeToFile(0);

		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.

		Real body_average_kinetic_energy = 0.0;
		Real body_maximum_kinetic_energy = 0.0;

		GlobalStaticVariables::physical_time_ = ite;
		while (ite < 25000)
		{
			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			body_inner.updateConfiguration();

			correction_matrix.exec();
			relaxation_inner.exec();

			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			ite++;

			if ((ite == 1) || (ite % 100 == 0))
			{
				update_body_kinetic_energy.exec();
				body_average_kinetic_energy = calculate_body_average_kinetic_energy.exec();
				body_maximum_kinetic_energy = calculate_body_maximum_kinetic_energy.exec();
				std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
				out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " << body_average_kinetic_energy << 
					"   " << body_maximum_kinetic_energy << "\n";
				write_body_to_plt.writeToFile(ite);
			}
		}
		ite++;

		write_body_to_plt.writeToFile(ite);
		write_particle_reload_files.writeToFile(0);
		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
}
