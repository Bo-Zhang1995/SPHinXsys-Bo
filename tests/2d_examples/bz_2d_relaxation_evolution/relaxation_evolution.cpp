/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case by testing the relaxation with evolution method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file
//----------------------------------------------------------------------
std::string input_body = "./input/TurbineBlade.dat";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 10.0;
Real DH = 10.0;
Real insert_circle_radius = 1.0;
Vec2d insert_circle_center(0.0, 0.0);
Real resolution_ref = 1 / 25.0;
Real BW = resolution_ref * 4;
BoundingBox system_domain_bounds(Vec2d(-BW - DL, -BW - DH), Vec2d(BW + DL, BW + DH));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
//class Circle : public ComplexShape
//{
//public:
//    explicit Circle(const std::string& shape_name) : ComplexShape(shape_name)
//    {
//        MultiPolygon multi_polygon;
//        multi_polygon.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
//        add<MultiPolygonShape>(multi_polygon);
//    }
//};

class InputBody : public ComplexShape
{
public:
    explicit InputBody(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon turbine_blade;
        turbine_blade.addAPolygonFromFile(input_body, ShapeBooleanOps::add);
        add<MultiPolygonShape>(turbine_blade);
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
    FluidBody body(sph_system, makeShared<InputBody>("Body"));
    body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    body.defineAdaptationRatios(0.8, 1.0);
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
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToPlt write_body_to_plt(io_environment, { &body });
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(io_environment, { &body });
        /** Calculate the correction matrix */
        InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> correction_matrix(body_inner);
        body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
        /** Relaxation method based on the constant pressure. */
        relax_dynamics::RelaxationStepInner relaxation_inner(body_inner, true);
        //relax_dynamics::RelaxationStepInnerImplicit relaxation_inner(body_inner, true);
        //relax_dynamics::RelaxationStepInner<CorrectionMatrixRelaxation> relaxation_inner(body_inner, true);
        //relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner(body_inner, true);
        /** Update the kinetic energy for the relaxation. */
        SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(body_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_insert_body_particles.exec(0.25);
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        write_body_to_plt.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        TickCount t1 = TickCount::now();
        int ite = 0;

        GlobalStaticVariables::physical_time_ = ite;
        while (ite < 100000)
        {
            correction_matrix.exec();
            relaxation_inner.exec();
            ite++;

            if (ite % 200 == 0)
            {
                update_body_kinetic_energy.exec();
                write_body_to_plt.writeToFile(ite);
                std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
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
};

